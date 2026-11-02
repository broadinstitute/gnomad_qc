"""Map which functions in each review area are called from other areas.

Used by the gnomad-pr-review skill to tell each reviewer agent who calls its
functions and which other areas it depends on.

Usage (run from the repo root of the PR being reviewed):

    python callmap.py FILE [FILE ...] [--split FILE:LINE[,LINE...]] [--callers GLOB ...]

- FILE: source files under review. Each file is one area.
- --split: break a large file into areas at the given line numbers, e.g.
  ``--split pkg/utils/constraint.py:825,1544`` gives areas [A] lines 1-824,
  [B] 825-1543, [C] 1544-end. Split at function boundaries. The file must be
  one of the FILEs.
- --callers: extra files (e.g. tests) that are scanned for calls but are not
  areas themselves. Default: every git-tracked ``.py`` file under the current
  directory (``**/*.py`` outside git). Pass ``--callers`` with no globs to
  scan only the FILEs.

Paths can be given in any spelling (``./x.py``, ``x.py``, absolute); they are
normalized relative to the current directory.

Output: one line per function or module-level name that is used from a
different area, with every caller as ``path:line(enclosing_function)``.
"""

import argparse
import ast
import collections
import glob
import os
import subprocess
import sys
from typing import Dict, List, Optional, Set, Tuple


def _norm(path: str) -> str:
    """Return a path normalized relative to the current directory."""
    return os.path.relpath(os.path.realpath(path))


def _defs(path: str) -> Dict[str, int]:
    """Return top-level function and assigned names in a file with their lines."""
    out = {}
    for node in ast.parse(open(path).read()).body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            out[node.name] = node.lineno
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    out.setdefault(t.id, node.lineno)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            out.setdefault(node.target.id, node.lineno)
    return out


def _module_stem(src: str, module: Optional[str], level: int) -> str:
    """
    Return the path (without ``.py``) that an import refers to.

    Absolute imports resolve from the current directory; relative imports
    (``level`` > 0) resolve from the importing file's package.
    """
    base = ""
    if level:
        base = os.path.dirname(src)
        for _ in range(level - 1):
            base = os.path.dirname(base)
    parts = module.split(".") if module else []
    return os.path.join(base, *parts) if (base or parts) else "."


def _dotted(node: ast.AST) -> Optional[str]:
    """Return ``a.b.c`` for a Name/Attribute chain, else None."""
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        inner = _dotted(node.value)
        return f"{inner}.{node.attr}" if inner else None
    return None


def _local_names(fn: ast.AST) -> Set[str]:
    """Return a function's parameters and assigned names, skipping nested scopes."""
    a = fn.args
    names = {
        x.arg for x in a.posonlyargs + a.args + a.kwonlyargs + [a.vararg, a.kwarg] if x
    }
    todo = list(ast.iter_child_nodes(fn))
    while todo:
        n = todo.pop()
        if isinstance(n, ast.Name) and not isinstance(n.ctx, ast.Load):
            names.add(n.id)
        elif not isinstance(
            n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Lambda)
        ):
            todo.extend(ast.iter_child_nodes(n))
    return names


def _default_callers() -> List[str]:
    """Return git-tracked .py files under cwd, or all .py files outside git."""
    try:
        out = subprocess.run(
            ["git", "ls-files", "--", "*.py"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        return out.split()
    except (OSError, subprocess.CalledProcessError):
        return glob.glob("**/*.py", recursive=True)


def main() -> None:
    """Parse arguments and print the cross-area call map."""
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("files", nargs="+")
    p.add_argument("--split", action="append", default=[])
    p.add_argument("--callers", nargs="*", default=None)
    a = p.parse_args()

    files = list(dict.fromkeys(_norm(f) for f in a.files))
    area_set = set(files)

    splits: Dict[str, List[int]] = {}
    for s in a.split:
        f, lines = s.rsplit(":", 1)
        f = _norm(f)
        if f not in area_set:
            p.error(f"--split file {f} is not one of the FILEs: {', '.join(files)}")
        splits[f] = sorted(int(x) for x in lines.split(","))

    # Label areas by basename, or by path when two FILEs share a basename.
    base_counts = collections.Counter(os.path.basename(f) for f in files)
    labels = {
        f: os.path.basename(f) if base_counts[os.path.basename(f)] == 1 else f
        for f in files
    }

    def _area(path: str, line: int) -> str:
        """Return the area label for a line in an area file."""
        if path not in splits:
            return labels[path]
        idx = sum(1 for b in splits[path] if line >= b)
        return f"{labels[path]}[{chr(ord('A') + idx)}]"

    # (defining file, name) -> line. Keyed by file so the same name defined in
    # two areas stays distinct; uses are resolved through imports.
    defs: Dict[Tuple[str, str], int] = {}
    for f in files:
        for name, line in _defs(f).items():
            defs[(f, name)] = line

    def _area_file(stem: str) -> Optional[str]:
        """Return the area file a module stem points to, if any."""
        for cand in (_norm(stem + ".py"), _norm(os.path.join(stem, "__init__.py"))):
            if cand in area_set:
                return cand
        return None

    callers = list(files)
    globs = (
        _default_callers()
        if a.callers is None
        else [c for g in a.callers for c in glob.glob(g, recursive=True)]
    )
    seen = set(callers)
    for c in globs:
        c = _norm(c)
        if c not in seen and c.endswith(".py"):
            seen.add(c)
            callers.append(c)

    uses = collections.defaultdict(list)
    for src in callers:
        try:
            tree = ast.parse(open(src).read())
        except (OSError, SyntaxError, UnicodeDecodeError) as e:
            print(f"# skipped {src}: {e}", file=sys.stderr)
            continue

        # Dotted module aliases (``import pkg.mod as alias``, ``from pkg import
        # mod``) that point at an area file, so ``alias.func()`` counts as a use
        # but ``args.func`` does not. Bare names count only if imported from an
        # area file (or used in the defining file itself), so a local variable
        # that shares a function's name in another file is not reported.
        mod_aliases: Dict[str, str] = {}  # local dotted name -> area file
        imported: Dict[str, Tuple[str, str]] = {}  # local name -> (file, name)
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for al in node.names:
                    f = _area_file(_module_stem(src, al.name, 0))
                    if f:
                        mod_aliases[al.asname or al.name] = f
            elif isinstance(node, ast.ImportFrom):
                stem = _module_stem(src, node.module, node.level)
                mod_file = _area_file(stem)
                for al in node.names:
                    if al.name == "*":
                        if mod_file:
                            print(
                                f"# {src}:{node.lineno}: star import from"
                                f" {mod_file}; its names are not tracked",
                                file=sys.stderr,
                            )
                        continue
                    sub_file = _area_file(os.path.join(stem, al.name))
                    if sub_file:
                        mod_aliases[al.asname or al.name] = sub_file
                    elif mod_file:
                        imported[al.asname or al.name] = (mod_file, al.name)

        class V(ast.NodeVisitor):
            def __init__(self) -> None:
                self.stack: List[Tuple[str, int]] = []
                self.local_names: List[Set[str]] = []  # per enclosing function

            def _scope(self, n: ast.AST) -> None:
                self.stack.append((n.name, n.lineno))
                self._visit_with_locals(n)
                self.stack.pop()

            def _visit_with_locals(self, n: ast.AST) -> None:
                is_cls = isinstance(n, ast.ClassDef)
                self.local_names.append(set() if is_cls else _local_names(n))
                self.generic_visit(n)
                self.local_names.pop()

            visit_FunctionDef = visit_AsyncFunctionDef = visit_ClassDef = _scope
            visit_Lambda = _visit_with_locals

            def _use(self, dfile: str, name: str, lineno: int) -> None:
                if (dfile, name) not in defs:
                    return
                dline = defs[(dfile, name)]
                # Attribute the use to the area of its enclosing top-level def.
                at = self.stack[0][1] if self.stack else lineno
                if src in area_set and _area(src, at) == _area(dfile, dline):
                    return
                where = "/".join(s[0] for s in self.stack) or "module"
                uses[(dfile, name)].append(f"{src}:{lineno}({where})")

            def visit_Name(self, n: ast.Name) -> None:
                # Skip assignments, and locals that share a function's name.
                if not isinstance(n.ctx, ast.Load) or any(
                    n.id in s for s in self.local_names
                ):
                    return
                if n.id in imported:
                    self._use(*imported[n.id], n.lineno)
                elif src in area_set:
                    self._use(src, n.id, n.lineno)

            def visit_Attribute(self, n: ast.Attribute) -> None:
                # Catches module-alias calls like ``constraint_res.get_x()``.
                mod = _dotted(n.value)
                if mod in mod_aliases:
                    self._use(mod_aliases[mod], n.attr, n.lineno)
                    return
                self.generic_visit(n)

        V().visit(tree)

    for (f, name), line in sorted(defs.items(), key=lambda kv: (kv[0][0], kv[1])):
        u = uses[(f, name)]
        if u:
            shown = ", ".join(u[:15]) + (
                f", ... (+{len(u) - 15})" if len(u) > 15 else ""
            )
            print(f"{_area(f, line)}:{line} {name}  <-  {shown}")


if __name__ == "__main__":
    main()
