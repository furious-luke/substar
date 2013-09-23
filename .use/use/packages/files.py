import use

class copy(use.Feature):
    def actions(self, sources, targets=[], options={}):
        return [use.actions.Copy(sources, targets)]

class run(use.Feature):
    def actions(self, sources, targets=[], options={}):
        return [use.actions.Command(sources.path, show_command=False, show_stdout=True)]

class Default(use.Version):
    version = 'default'
    features = [copy, run]

class files(use.Package):
    default_target_node = use.File
    default_builder = use.Builder
    versions = [Default]
