# Copyright © 2021 Jonathan Grot
#
# This file is part of grotlib. The latest version of grotlib is
# available at https://github.com/jgrot/grotlib
#
# Grotlib is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Grotlib is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Grotlib.  If not, see <https://www.gnu.org/licenses/>.
#

import collections
import six

def isIterableNonString( item ) :
    """From https://stackoverflow.com/questions/1055360/how-to-tell-a-variable-is-iterable-but-not-a-string

    See Answer from sorin.

    TODO: find new home.
    """

    return (
        isinstance(item, collections.Iterable) 
        and not isinstance(item, six.string_types)
    )
