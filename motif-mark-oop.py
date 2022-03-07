'''
Beagan Nguy
OOP Motif Mark Assignment
Bi 625

Usage: python motif-mark-oop.py -f <fasta file> -m <text file containing a list of motifs (1 motif per line)>
Output: Visualization of motifs on sequences as PNG file

Note:
- Motifs are translucent to show overlap
- Program can scale to infinite motifs and infinite sequences as inputs
'''
import argparse
import sys
import regex
import re
import pprint
import cairo
import random

# Fasta file handler. Return header and sequence
class FastAreader():
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta(self):
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split())
						
        yield header, sequence

# Command line handler
class CommandLine():
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Input fasta file and motif file'
        )
 
        self.parser.add_argument(
            '-f', '-fasta', type=str, help='Specify fasta file', required=False
        )

        self.parser.add_argument(
            '-m', '-motif_file', help='Specify the list of motif file', required=False
        )

        self.args = self.parser.parse_args()

# Distinguishes Exon and Intron of a sequence by its capitalization.
# low case = intron, upper case = exon
class ExonIntron():
    def __init__(self, header, sequence, count):
        self.header = header
        self.sequence = sequence
        self.key = header.split(" ")[0] + "-{}".format(count)
        self.count = count
        self.exon_intron_list_strip = self.findExonIntron()
    
    def findExonIntron(self):
        exon_intron_list = re.findall('[a-z]*|[A-Z]*', self.sequence.strip(" "))
        exon_intron_list_strip = [i for i in exon_intron_list if i]
        return(exon_intron_list_strip)

    def returnKey(self):
        return(self.key)

# Ambigous nucleotide motif handler
# Dictionary and regex was implemented to search for motif coordinates
class MotifDegenerate():
    def __init__(self, motif_list):
        self.motif_list = motif_list
        self.degenerate_conversion_upper = {"W":"[A|T]", "S":"[C|G]", "M":"[A|C]", "K":"[G|T]", "R":"[A|G]", "Y":"[C|T]",
                                      "B":"[C|G|T]", "D":"[A|G|T]", "H":"[A|C|T]", "V":"[A|C|G]",
                                      "N":"[A|C|G|T]",
                                      "A":"[A|U]", "C":"[C]", "U":"[U]", "G":"[G]", "T":"[T]"
                                    }
        self.degenerate_conversion_lower = {"w":"[a|t]", "s":"[c|g]", "m":"[a|c]", "k":"[g|t]", "r":"[a|g]", "y":"[c|t]",
                                      "b":"[c|g|t]", "d":"[a|g|t]", "h":"[a|c|t]", "v":"[a|c|g]",
                                      "n":"[a|c|g|t]",
                                      "a":"[a|u]", "c":"[c]", "u":"[u]", "g":"[g]", "t":"[t]"
                                    }
        self.motif_permutation_list = self.buildRegexMotif()

    def buildRegexMotif(self):
        print(self.motif_list)
        a_motif_list = list()
        for motif in self.motif_list:
            a_motif = str()
            for nucleotide in motif:   
                if nucleotide.islower():
                    motif_regex = self.degenerate_conversion_lower[nucleotide] 
                    a_motif += motif_regex
                else:
                    motif_regex = self.degenerate_conversion_upper[nucleotide] 
                    a_motif += motif_regex
            a_motif_list.append(a_motif)
        return (a_motif_list)
    
    def returnMotif(self):
        return(self.motif_permutation_list)

# Use pycairo to draw motifs on sequences
class Draw():
    def __init__(self, exon_intron, matches, motif_list, header_list, matches_dict, fasta_name):
        self.fasta_name = fasta_name
        self.exon_intron = exon_intron
        image_scaler = (len(exon_intron) // 10) + 1
        self.image_number = 1000*image_scaler
        self.matches = matches
        self.motif_list = motif_list
        self.header = header_list
        self.matches_dict = matches_dict
        self.draw_track = cairo.ImageSurface(cairo.FORMAT_RGB24, self.image_number + 200, self.image_number + 200)
        self.color_motif_dict = dict()
        self.colorMotifAssign()
        self.drawBackground()
        self.drawCap()
        self.drawExonIntron()
        self.drawMotif()
        self.drawLegend()
        self.convertToPNG()
    
    # Randomly assign colors to motifs
    def colorMotifAssign(self):
        for head, motif in self.matches_dict.items():
            for key, value in motif.items():
                r = random.randint(0, 100)/100
                g = random.randint(0, 100)/100
                b = random.randint(0, 100)/100
                self.color_motif_dict[key] = [r, g, b]
            break
    
    # Draw the background
    def drawBackground(self):
        context = cairo.Context(self.draw_track)
        context.rectangle(0, 0, self.image_number + 200, self.image_number + 200)
        context.set_source_rgb(1, 1, 1)
        context.fill()

    # Draw end caps
    def drawCap(self):
        context = cairo.Context(self.draw_track)
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(16)

        for count in range(1, len(self.exon_intron) + 1):
            context.set_line_cap(cairo.LINE_CAP_BUTT)
            context.move_to(30, ((self.image_number/len(self.exon_intron))*count) - 120)
            context.line_to(35, ((self.image_number/len(self.exon_intron))*count) - 120)
            context.stroke()
    
    # Draw exons and introns
    def drawExonIntron(self):
        context = cairo.Context(self.draw_track)

        for count, seq_list in enumerate(self.exon_intron):
            startset_off = 35
            count += 1
            for seq in seq_list:
                # Lower = introns
                if seq[0].islower():
                    #print("[{}]".format(seq))
                    context.set_source_rgba(0.5, 0.5, 0.60, 0.9)

                    # Setting of line width
                    context.set_line_width(6)

                    # Move the context to x,y position
                    context.move_to(startset_off, ((self.image_number/len(self.exon_intron))*count) - 120)

                    # Creating a line
                    #print(len(seq) + startset_off)
                    context.line_to(len(seq) + startset_off, ((self.image_number/len(self.exon_intron))*count) - 120)

                    # Stroke out the color and width property
                    context.stroke()

                    # Offset start
                    startset_off += len(seq)
                
                # Upper = Exon
                if seq[0].isupper():
                    #print("[{}]".format(seq))
                    context.set_source_rgba(0, 0, 0.60, 0.7)

                    # Setting of line width
                    context.set_line_width(18)

                    # Move the context to x,y position
                    context.move_to(startset_off, ((self.image_number/len(self.exon_intron))*count) - 120)

                    # Creating a line
                    #print(len(seq) + startset_off)
                    context.line_to(len(seq) + startset_off, ((self.image_number/len(self.exon_intron))*count) - 120)

                    # Stroke out the color and width property
                    context.stroke()

                    # Offset start
                    startset_off += len(seq)
            
            # Draw 3' cap
            context.set_source_rgb(0, 0, 0)
            context.set_line_width(16)

            context.set_line_cap(cairo.LINE_CAP_BUTT)
            context.move_to(startset_off, ((self.image_number/len(self.exon_intron))*count) - 120)
            context.line_to(startset_off+5, ((self.image_number/len(self.exon_intron))*count) - 120)
            context.stroke()

            # Draw label
            context.set_source_rgb(0.1, 0.1, 0.1)
            context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            context.set_font_size(14)

            context.move_to(30, ((self.image_number/len(self.exon_intron))*count) - 170)
            head = self.header[count - 1]
            context.show_text(head)

    # Draw motifs  
    def drawMotif(self):
        context = cairo.Context(self.draw_track)
        startset_off = 35
        head_count = 0
        for head, motif_keys in self.matches_dict.items():
            head_count += 1
            for motif, coord_set in motif_keys.items():
                #print(coord_set)
                for coord in coord_set:
                    if coord:
                        color = self.color_motif_dict[motif]

                        context.set_source_rgba(color[0], color[1], color[2], 0.80)

                        # Setting of line width
                        context.set_line_width(24)

                        # Move the context to x,y position
                        context.move_to(startset_off + int(coord[0]), ((self.image_number/len(self.exon_intron))*head_count) - 120)

                        # Creating a line
                        #print(len(seq) + startset_off)
                        context.line_to(startset_off + int(coord[1]), ((self.image_number/len(self.exon_intron))*head_count) - 120)

                        # Stroke out the color and width property
                        context.stroke()
    
    # Draw legend
    def drawLegend(self):
        space = 0
        context = cairo.Context(self.draw_track)
        for motif, color in self.color_motif_dict.items():
            context.set_source_rgba(color[0], color[1], color[2], 0.70)
            context.set_line_width(10)
            context.move_to(self.image_number, 100 + space)
            context.line_to(self.image_number + 150, 100 + space)
            space += 10
            context.stroke()

            # Make blank line
            context.set_source_rgb(1, 1, 1)
            context.set_line_width(10)
            context.move_to(self.image_number, 100 + space + 10)
            context.line_to(self.image_number + 150, 100 + space + 10)
            space += 20
            context.stroke()

             # Draw label
            context.set_source_rgb(0, 0, 0)
            context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            context.set_font_size(14)
            context.move_to(self.image_number, 100 + space - 40)
            context.show_text(motif)
    
    # SVG to PNG
    def convertToPNG(self):
        self.draw_track.write_to_png('{}.png'.format(self.fasta_name[:-6]))

# Main Function 
def main():
    myCommandLine = CommandLine()
    reader = FastAreader(myCommandLine.args.f)
    count = 0 # Count to be used to count the number of figures
    fh = open(myCommandLine.args.m, "r")
    og_motif_list = [motif.strip("\n") for motif in fh]
    motif_list = MotifDegenerate(og_motif_list).returnMotif()
    matches_dict = dict()
    matches = list() # Will contain list of motif coordinate matches
    
    # Initiate a bunch of list in a list. To be used to hold list of motifs and their coordinate matches
    for i in range(0, len(motif_list)):
        matches.append(list())
    
    exon_intron_list = list()
    header_list = list()

    # Iterate through each header and sequence pair in a fasta file
    for head, sequence in reader.readFasta():
        header_list.append(head)
        temp_matches = list()
        count += 1
        myExon = ExonIntron(head, sequence, count)
        coord = myExon.exon_intron_list_strip
        exon_intron_list.append(coord)
        matches_dict[head] = dict()
        #print("--{}--".format(motif_list))
        # Get coordinates for motifs on sequences. Overlapped set to true to find overlapping motifs
        for count, motif in enumerate(motif_list):
            temp_matches = ([(m.start(0), m.end(0)) for m in regex.finditer(motif, sequence, overlapped=True)])
            matches[count].append(temp_matches)
            matches_dict[head][og_motif_list[count]] = temp_matches
    #pprint.pprint(matches[0])
    #pprint.pprint(exon_intron_list[0])
    #pprint.pprint(matches_dict)
    Draw(exon_intron_list, matches, motif_list, header_list, matches_dict, myCommandLine.args.f)
if __name__ == "__main__":
	main()