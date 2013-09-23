from subprocess import Popen, PIPE

class Tool(object ):

    def __init__(self, cmd, args):
        self.cmd = cmd
        self._args

    def run(self):
        proc = Popen(self._cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate()
        return (proc.returncode, stdout, stderr)

    def _setup_cmd_line(self):
        return '%s ' + self._args + self.name self._cmd
