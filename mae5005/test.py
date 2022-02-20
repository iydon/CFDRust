import json
import pathlib as p
import re
import subprocess as s
import typing as t

try:
    import colorama

    colorama.init()
except ImportError:
    colorama = None


Path = t.Union[str, p.Path]


class Fortran:
    def __init__(self, path: Path) -> None:
        self._path = p.Path(path)
        self._is_compiled = False

    def compile(self, *args: str) -> 'Fortran':
        if not self._is_compiled:
            self._is_compiled = self._compile(args)
        return self

    def test(self) -> bool:
        self.compile()
        flag = True
        for test in self._tests():
            stdin, stdout = test.get('stdin', None), test.get('stdout', None)
            if self._run(stdin, stdout):
                print(self._color('[pass]', 'green'), end=' ')
            else:
                print(self._color('[fail]', 'red'), end=' ')
                flag = False
            print(f'stdin={repr(stdin)}, stdout={repr(stdout)}')
        return flag

    def _compile(self, args: t.Tuple[str, ...]) -> bool:
        command = self._command()
        args = [
            command.get('fortran', 'gfortran'), *command.get('options', []),
            str(self._path), '-o', f'{self._path.stem}.out', *args,
        ]
        try:
            return s.run(args).returncode == 0
        except Exception as e:
            print(f'{e}: {" ".join(args)}')
            return False

    def _run(self, stdin: t.Optional[str] = None, stdout: t.Optional[str] = None) -> bool:
        with s.Popen(f'./{self._path.stem}.out', stdin=s.PIPE, stdout=s.PIPE) as process:
            input = None if stdin is None else stdin.encode()
            out, err = process.communicate(input)
            if err is not None:
                print(err.decode())
                return False
            elif stdout is not None:
                return out.decode().strip() == stdout.strip()
            return True

    def _tests(self) -> t.List[t.Dict[str, str]]:
        try:
            return json.loads('\n'.join(self._label_extractor('#[test]')))
        except json.decoder.JSONDecodeError:
            return []

    def _command(self) -> t.Dict[str, t.Any]:
        try:
            return json.loads('\n'.join(self._label_extractor('#[compile]')))
        except json.decoder.JSONDecodeError:
            return {}

    def _label_extractor(self, annotation: str) -> t.Iterator[str]:
        is_test = False
        for line in self._path.read_text().splitlines():
            if is_test:
                if line.startswith('!'):
                    yield line[1:]
                else:
                    is_test = False
            elif line == f'! {annotation}':
                is_test = True

    def _color(self, string: str, color: str) -> str:
        if colorama is None:
            return string
        return getattr(colorama.Fore, color.upper()) + string + colorama.Style.RESET_ALL


if __name__ == '__main__':
    import os
    import sys

    fails = []
    columns = os.get_terminal_size().columns
    for arg in sys.argv[1:]:
        print(f'{arg:-^{columns}}')
        if not Fortran(arg).compile('-O3').test():
            fails.append(arg)
    print(f'{"FailedCode":=^{columns}}')
    print('[', ', '.join(fails), ']')
