import os, tempfile, shutil, unittest
from ..Graph import Graph

class Rule(object):

    def productions(self, node):
        return [((0,), 7), ((1,), 8), ((2,), 9)]

class TestGraph(unittest.TestCase):

    def test_do_expand_node(self):
        g = self.graph._graph
        for ii in range(7):
            g.add_node(ii)
        g.add_edge(0, 3)
        g.add_edge(1, 3)
        g.add_edge(2, 3)
        g.add_edge(3, 4)
        g.add_edge(3, 5)
        g.add_edge(3, 6)
        self.graph._do_expand_node(3, Rule())
        self.assertSetEqual(set(g.successors(0)), set([7]))
        self.assertSetEqual(set(g.successors(1)), set([8]))
        self.assertSetEqual(set(g.successors(2)), set([9]))
        self.assertSetEqual(set(g.successors(7)), set([4, 5, 6]))
        self.assertSetEqual(set(g.successors(8)), set([4, 5, 6]))
        self.assertSetEqual(set(g.successors(9)), set([4, 5, 6]))

    def setUp(self):
        self.graph = Graph()

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
