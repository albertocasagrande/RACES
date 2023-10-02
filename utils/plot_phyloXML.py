#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Alberto Casagrande <alberto.casagrande@uniud.it>
#
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sys
import re
import os.path
import argparse
from ete3 import Phyloxml, NodeStyle


def from_branchcolor_to_color(branchcolor):
    return '#' + (hex(branchcolor.red)[2:].rjust(2, '0')
                  + hex(branchcolor.green)[2:].rjust(2, '0')
                  + hex(branchcolor.blue)[2:].rjust(2, '0'))


def image_file_type(filename, pat=re.compile(r"^.*\.png|.*\.svg$")):
    if os.path.exists(filename):
        raise argparse.ArgumentTypeError(f'"{filename}" already exists')

    if not pat.match(filename):
        raise argparse.ArgumentTypeError('the image file extension '
                                         + 'must be either "svg" or'
                                         + ' "png"')
    return filename


if __name__ == '__main__':
    description = 'A simple program to plot a phyloXML file'

    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description=description)
    parser.add_argument('XMLfile', type=argparse.FileType('r'),
                        help='the input phyXML filename')
    parser.add_argument('-o', '--output-image', dest='image_filename',
                        help=('an optional output image filename'
                              + ' (png or svg)'),
                        type=image_file_type)

    args = parser.parse_args()

    project = Phyloxml()
    project.build_from_file(args.XMLfile.name)

    trees = [tree for tree in project.get_phylogeny()]

    if args.image_filename is not None:
        suffix = args.image_filename[-3:]
        prefix = args.image_filename[:-4]

    i = 0
    for tree in trees:
        for node in tree.traverse():
            branchcolor = node.phyloxml_clade.get_color()

            nstyle = NodeStyle()
            if branchcolor is not None:
                color = from_branchcolor_to_color(branchcolor)
                nstyle["fgcolor"] = color
            nstyle["size"] = 5
            node.set_style(nstyle)

        tree.ladderize()

        tree.show()

        if args.image_filename is not None:
            if len(trees) > 1:
                tree.render(f'{prefix}_{i}.{suffix}')
                i += 1
            else:
                tree.render(args.image_filename)
