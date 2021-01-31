# Copyright © 2021 Jonathan Grot

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
