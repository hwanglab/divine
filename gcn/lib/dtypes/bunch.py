"""
A Bunch is a generic class to which one can dynamically add attributes.
It supports setting and accessing of attributes using both dictionary-style
access as well as attribute-style access. It also supports the keys,
values and items methods of a dictionary.

Usage

>>> stuff = Bunch(grapes=20, bananas=10)
>>> stuff.grapes
20
>>> stuff.grapes = 5
>>> stuff.grapes
5
>>> stuff.bananas
10
>>> stuff['bananas'] = 21
>>> stuff.bananas
21
>>> stuff.apples
False
>>> stuff.pineapples = 9
>>> stuff['pineapples'] + stuff['bananas']
30
>>> alist = list(stuff.keys())
>>> alist.sort()
>>> print(alist)
['bananas', 'grapes', 'pineapples']
"""


class Bunch(dict):

    def __getattr__(self, v):
        return self.get(v, False)

    def __setattr__(self, v, item):
        self[v] = item

    def __delattr__(self, v):
        try:
            self.pop(v)
        except KeyError:
            pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
