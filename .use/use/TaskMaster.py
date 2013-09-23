class TaskMaster(object):

    def __init__(self, graph):
        self.graph = graph

    def __call__(self):
        targets = graph.targets()
        for node in targets:
            self._build_node(node)

    def _build_node(self, node):
        if self.graph.is_node_built(node):
            return
        parents = graph.parents(node)
        for parent in parents:
            self._build_node(parent)
        self._do_node_build(node)

    def _do_node_build(self, node):
        pass
