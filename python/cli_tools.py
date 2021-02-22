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
import re
import six
import sys
import tabulate

# Grotlib imports
import type_tools

def addVerboseOptionTo( subparser ) :
    """Adds a verbose (dest=debug) option to an argparse parser

    """
    subparser.add_argument('-v', dest='debug', default=False, action='store_true',
                           help='Print verbose info for debugging purposes')


def cli_closure( closure, dump_traceback = False ) :
    """Convenience wrapper for running CLI commands.

    This routine nicely traps Exceptions, formats the exception
    message, writes to STDERR and calls exit( 1 )

    Args:

       closure (lambda): The function and environment to be executed.
       dump_traceback (bool): If True, dumps the whole traceback as would normally be done.  Does not try to trap exceptions.  Useful while developing the closure.

    Usage:

       cli_closure( lambda : my_cli_command( args ) )

    TODO:find new home?

    """

    if dump_traceback :

        closure( )

    else :

        errmsg = None

        try :
            
            closure( )
            
        except Exception as e :

            # One would think work would not be needed to extract the
            # final error string from an exception, but I can't find any
            # examples showing how to do that.  So here's my version.
            
            errtyp = str(type(e))
            
            # Extract the Exception class name
            m = re.match( "\<class '(.*)'\>", errtyp)
            fullclassname = m.group(1)
            
            # Take the last dotted field
            classname = fullclassname.split('.')[-1]

            errmsg = "%s: %s" % ( classname, str(e) )

        if errmsg is not None :

            sys.stderr.write( "***Error: %s\n" % errmsg )

            exit( 1 )


class RowWriter :
    """Writes first N columns with fixed widths and the N+1'th column any width

    Also has a divider line writer.  See the __main__ section for an
    example.  Also writes a tabulate table.

    """
    
    def __init__( self, column_widths = [20], stream = sys.stdout,
                  justification = None, indent = 0, sep = ": " ) :

        if any( map( lambda w: w<=0, column_widths ) ) :
            raise Exception( "Bad width in column_widths" )

        self.ncols = len( column_widths )
        self.column_widths = column_widths
        self.indent = indent
        self.sep = sep
        self.stream = stream

        # This grows as the maximum length string in write grows.
        # Limit is DIVMAX.
        self.divider_length = 0
        self.DIVMAX = 80

        if justification is not None :

            if len( justification ) != self.ncols :
                raise Exception( "Justification must contain 'l' 'r' or 'c' for each column" )
            for c in justification :
                if c not in ["l","r","c"] :
                    raise Exception( "Bad justification character: %s.  Must be 'l' 'r' or 'c'" % c )

            self.justification = justification

        else :

            self.justification = "r" * (self.ncols)


    def divider( self, char = "-" ) :
        """Write a divider line.

        The divider line grows with the longest line printed up to the
        limit self.DIVMAX.

        """

        self.stream.write( " "*self.indent )

        if self.divider_length == 0 :

            self.stream.write( char*self.DIVMAX )

        else :

            self.stream.write( char*self.divider_length )

        self.stream.write( "\n" )
        
    def print( self, s ) :
        """Writes the string to self.stream and adds a newline
        """

        lines = s.split( '\n' )
        
        for line in lines :

            self.divider_length = min( max( self.divider_length, len(line) ), self.DIVMAX )
            self.stream.write( " "*self.indent )
            self.stream.write( line )
            self.stream.write( '\n' )


    def tabulate( self, headers, data ) :
        """Calls tabulate.tabulate on headers and data"""

        self.print( tabulate.tabulate( data, headers ) )


    def write( self, *items, list_in_col = False ) :
        """Write items in a table of sorts.

        The last column is variable width and left justified.  There
        must be one more item than column width specifications when
        creating a RowWriter object.

        if list_in_col is true, if the last item is a list the list
        will be written on multple lines, one line per list element.

        """

        if len(items) != self.ncols + 1 :
            raise Exception('Expected an item count of one plus the number of fixed width columns')

        self.stream.write( ' '*self.indent )

        s4 = ""
        
        for iitem, item in enumerate( items ) :

            if iitem == self.ncols :
                
                if type_tools.isIterableNonString( item ) and list_in_col :

                    s3 = repr( item[0] )
                    keep_writing = item[1:]

                else :

                    if not isinstance( item, six.string_types ) :
                        s3 = repr( item )
                    else :
                        s3 = item

                    keep_writing = []

            else :

                if not isinstance( item, six.string_types ) :
                    # Not a string, convert to string
                    s1 = repr( item )
                else :
                    s1 = item

                w = self.column_widths[iitem]

                if len(s1) > w :

                    s2 = s1[:(w-1)] + "*"

                else :

                    s2 = s1

                j = self.justification[iitem]

                if j == "r" :
                    fmt = "{:>%i}" % w
                    s3 = fmt.format( s2 )
                elif j == "l" :
                    fmt = "{:<%i}" % w
                    s3 = fmt.format( s2 )
                else :
                    # center
                    fmt = "{:^%i}" % w
                    s3 = fmt.format( s2 )

                s3 += self.sep
                
            s4 += s3

        self.divider_length = min( max( self.divider_length, len(s4) ), self.DIVMAX )
        
        s4 += "\n"        

        self.stream.write(s4)

        if len( keep_writing ) > 0 :
            
            offset_width = sum( self.column_widths ) + len( self.sep )*self.ncols

            for item in keep_writing :
                # Note: s4 does not include the indent
                s4 = " "*offset_width + repr(item)
                self.divider_length = min( max( self.divider_length, len(s4) ), self.DIVMAX )
                s5 = " "*self.indent + s4 + '\n'
                self.stream.write( s5 )

if __name__ == "__main__" :

    print("Demonstrating RowWriter")

    rw = RowWriter( [10, 6], sys.stdout, indent=4 )
    
    rw.divider( )
    rw.write( "action", "return", 42 )
    rw.divider( )
    rw.write( "var", "a", .24 )
    rw.divider( )
    rw.write( "other thing", "?", "xxxxx")
    rw.divider( )
    rw.write( "list of things", "the list", [1,2,3,4,5,"a long string.  A very very very very long string"], list_in_col=True)
    rw.divider( )
