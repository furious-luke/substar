import os, subprocess, shlex, shutil

class Installer(object):

    def __init__(self):
        pass

    def download_package(self, url, path):

        # Remove existing file.
        if os.path.exists(path):
            os.remove(path)

        # Attempt to download the package.
        try:
            import urllib
            urllib.urlretrieve(url, path)
            return True
        except:
            return False

    def unpack_package(self, src_path, dst_path):

        # Remove the destination entirely.
        if os.path.exists(dst_path):
            shutil.rmtree(dst_path)

        if os.path.splitext(src_path)[1] == '.zip':
            try:
                import zipfile
                zf = zipfile.ZipFile(src_path)
                zf.extractall(dst_path)
                zf.close()
            except:
                shutil.rmtree(dst_path, True)
                return False
        else:
            try:
                import tarfile
                tf = tarfile.open(src_path)
                tf.extractall(dst_path)
            except:
                shutil.rmtree(dst_path, True)
                return False

    def build(self, pkg, dst_path):

        # Remove any existing file used to indicate successful builds.
        if os.path.exists('use_build_success'):
            os.remove('use_build_success')

        # Remove the installation directory.
        if os.path.exists(dst_path):
            shutil.rmtree(dst_path)

        # Hunt down the correct build handler.
        handler = pkg.build_handler()
        if handler is None:
            return False

        # Make a file to log stdout from the commands.
        with open('stdout.log', 'w') as stdout_log:

            # Process each command in turn.
            for cmd in handler:

                # It's possible to have a tuple, indicating a function and arguments.
                if isinstance(cmd, tuple):
                    ctx.Log("Command is a Python function\n")
                    func = cmd[0]
                    args = cmd[1:]

                    # Perform substitutions.
                    args = [ctx.env.subst(a.replace('${PREFIX}', dst_path)) for a in args]

                    # Call the function.
                    func(*args)

                else:

                    # If the first character in a command is an "!", then it means we allow
                    # errors from this command.
                    allow_errors = False
                    if cmd[0] == '!':
                        allow_errors = True
                        cmd = cmd[1:]

                    # Perform substitutions.
                    cmd = cmd.format(prefix=dst_path)

                    try:
                        subprocess.check_call(shlex.split(cmd), stdout=stdout_log, stderr=subprocess.STDOUT)
                    except:
                        if not allow_errors:
                            return False

        # If it all seemed to work, write a dummy file to indicate this package has been built.
        success = open('use_build_success', 'w')
        success.write(' ')
        success.close()

        return True
