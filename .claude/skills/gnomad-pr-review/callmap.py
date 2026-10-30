"""Map which functions in each review area are called from other areas.

Used by the gnomad-pr-review skill to tell each reviewer agent who calls its
functions and which other areas it depends on.

Usage (run from the repo root of the PR being reviewed):

    python callmap.py FILE [FILE ...] [--split FILE:LINE[,LINE...]] [--callers GLOB ...]

- FILE: source files under review. Each file is one area.
- --split: break a large file into areas at the given line numbers, e.g.
  ``--split pkg/utils/constraint.py:825,1544`` gives areas [A] lines 1-824,
  [B] 825-1543, [C] 1544-end. Split at function boundaries.
- --callers: extra files (e.g. tests) that are scanned for calls but are not
  areas themselves. Default: ``tests/**/*.py``.

Output: one line per function or module-level name that is used from a
different area, with every caller as ``file:line(enclosing_function)``.
"""

import argparse
import ast
import collections
import glob
import os
import sys
from typing import Dict, List, Tuple


def _area(path: str, line: int, splits: Dict[str, List[int]]) -> str:
    """Return the area label for a line in a file."""
    name = os.path.basename(path)
    if path not in splits:
        return name
    bounds = splits[path]
    idx = sum(1 for b in bounds if line >= b)
    return f"{name}[{chr(ord('A') + idx)}]"


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


def main() -> None:
    """Parse arguments and print the cross-area call map."""
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("files", nargs="+")
    p.add_argument("--split", action="append", default=[])
    p.add_argument("--callers", nargs="*", default=["tests/**/*.py"])
    a = p.parse_args()

    splits: Dict[str, List[int]] = {}
    for s in a.split:
        f, lines = s.rsplit(":", 1)
        splits[f] = sorted(int(x) for x in lines.split(","))

    # name -> (defining file, line). Names defined in several files are skipped,
    # since an unqualified use can't be attributed.
    defs: Dict[str, Tuple[str, int]] = {}
    dupes = set()
    for f in a.files:
        for name, line in _defs(f).items():
            if name in defs:
                dupes.add(name)
            defs[name] = (f, line)
    for d in dupes:
        defs.pop(d, None)

    callers = list(a.files)
    for g in a.callers:
        callers += [c for c in glob.glob(g, recursive=True) if c not in callers]

    uses = collections.defaultdict(list)
    for src in callers:
        try:
            tree = ast.parse(open(src).read())
        except (SyntaxError, UnicodeDecodeError) as e:
            print(f"# skipped {src}: {e}", file=sys.stderr)
            continue

        # Module aliases (``import pkg.mod as alias``) that point at an area file,
        # so ``alias.func()`` counts as a use but ``args.func`` does not.
        # Bare names count only if imported from the defining area file (or used in
        # the defining file itself), so a local variable that shares a function's
        # name is not reported.
        area_paths = {os.path.normpath(f) for f in a.files}
        mod_aliases = set()
        imported = {}  # local alias -> original name
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for al in node.names:
                    path = os.path.normpath(al.name.replace(".", "/") + ".py")
                    if path in area_paths:
                        mod_aliases.add(al.asname or al.name.split(".")[-1])
            elif isinstance(node, ast.ImportFrom) and node.module:
                mod_path = os.path.normpath(node.module.replace(".", "/") + ".py")
                for al in node.names:
                    sub_path = os.path.normpath(
                        f"{node.module.replace('.', '/')}/{al.name}.py"
                    )
                    if sub_path in area_paths:
                        mod_aliases.add(al.asname or al.name)
                    elif mod_path in area_paths:
                        imported[al.asname or al.name] = al.name

        class V(ast.NodeVisitor):
            def __init__(self) -> None:
                self.stack: List[Tuple[str, int]] = []

            def _scope(self, n: ast.AST) -> None:
                self.stack.append((n.name, n.lineno))
                self.generic_visit(n)
                self.stack.pop()

            visit_FunctionDef = visit_AsyncFunctionDef = visit_ClassDef = _scope

            def _use(self, name: str, lineno: int) -> None:
                if name not in defs:
                    return
                dfile, dline = defs[name]
                # Attribute the use to the area of its enclosing top-level def.
                at = self.stack[0][1] if self.stack else lineno
                if src in a.files and _area(src, at, splits) == _area(
                    dfile, dline, splits
                ):
                    return
                where = "/".join(s[0] for s in self.stack) or "module"
                uses[name].append(f"{os.path.basename(src)}:{lineno}({where})")

            def visit_Name(self, n: ast.Name) -> None:
                if n.id in imported:
                    self._use(imported[n.id], n.lineno)
                elif n.id in defs and os.path.normpath(defs[n.id][0]) == (
                    os.path.normpath(src)
                ):
                    self._use(n.id, n.lineno)

            def visit_Attribute(self, n: ast.Attribute) -> None:
                # Catches module-alias calls like ``constraint_res.get_x()``.
                if isinstance(n.value, ast.Name) and n.value.id in mod_aliases:
                    self._use(n.attr, n.lineno)
                self.generic_visit(n)

        V().visit(tree)

    for name, (f, line) in sorted(
        defs.items(), key=lambda kv: (kv[1][0], kv[1][1])
    ):
        if uses[name]:
            u = uses[name]
            shown = ", ".join(u[:15]) + (f", ... (+{len(u) - 15})" if len(u) > 15 else "")
            print(f"{_area(f, line, splits)}:{line} {name}  <-  {shown}")


if __name__ == "__main__":
    main()
