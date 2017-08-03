import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
parser.add_argument('-ws', type=int, metavar='window_size', required=False, default='50000', help='Size of windows in bp')
args = parser.parse_args()

win_start = 0
win_end = args.ws
oldscaff = "scaffold_1"
num_CDS_bases = 0
out1 = open(args.o + "GeneDensity_WS" + str(args.ws / 1000) + "k.txt", 'w')
out1.write("scaff\tstart\tend\tGD\n")
with open(args.i, 'r') as gff:
    for i, line in enumerate(gff):
        # print(line)
        line = line.strip("\n").split("\t")
        if line[0].split("_")[0] == "scaffold":
            scaff, version, entry_type, start, stop, a1, a2, a3, info = line
            start = int(start)
            stop = int(stop)
            if entry_type == "CDS":

                if all(k > win_start for k in [start, stop]) and all(k <= win_end for k in [start, stop]) and scaff == oldscaff:  # CDS is entirely within window bounds
                    num_CDS_bases += stop - start + 1
                if start > win_start and start < win_end and stop > win_end and scaff == oldscaff:  # CDS is overhanging end of window
                    # Move onto next window, but parse overhang appropriately
                    num_CDS_bases += win_end - start + 1


                    GD = float(num_CDS_bases) / float(args.ws)

                    out1.write(scaff + "\t" +
                               str(win_start) + "\t" +
                               str(win_end) + "\t" +
                               str(GD) + "\n")

                    num_CDS_bases = stop - win_end
                    win_start = win_end
                    win_end += args.ws

                elif start > win_end and scaff == oldscaff:
                    GD = float(num_CDS_bases) / float(args.ws)

                    out1.write(scaff + "\t" +
                               str(win_start) + "\t" +
                               str(win_end) + "\t" +
                               str(GD) + "\n")

                    while start > win_end:
                        win_end += args.ws
                        win_start = win_end - args.ws
                        if start > win_end:
                            out1.write(scaff + "\t" +
                                       str(win_start) + "\t" +
                                       str(win_end) + "\t" +
                                       str(0.0) + "\n")

                    num_CDS_bases = 0

                    if all(k > win_start for k in [start, stop]) and all(k <= win_end for k in [start, stop]) and scaff == oldscaff:  # CDS is entirely within window bounds
                        num_CDS_bases += stop - start + 1
                    if start > win_start and start < win_end and stop > win_end and scaff == oldscaff:  # CDS is overhanging end of window
                        # Move onto next window, but parse overhang appropriately
                        num_CDS_bases += win_end - start + 1

                        GD = float(num_CDS_bases) / float(args.ws)

                        out1.write(scaff + "\t" +
                                   str(win_start) + "\t" +
                                   str(win_end) + "\t" +
                                   str(GD) + "\n")

                        num_CDS_bases = stop - win_end
                        win_start = win_end
                        win_end += args.ws
                    

                elif scaff != oldscaff:
                    oldscaff = scaff

                    win_start = 0
                    win_end = args.ws
                    num_CDS_bases = 0

                    if all(k > win_start for k in [start, stop]) and all(k <= win_end for k in [start, stop]) and scaff == oldscaff:  # CDS is entirely within window bounds
                        num_CDS_bases += stop - start + 1
                    if start > win_start and start < win_end and stop > win_end and scaff == oldscaff:  # CDS is overhanging end of window
                        # Move onto next window, but parse overhang appropriately
                        
                        num_CDS_bases += win_end - start + 1

                        GD = float(num_CDS_bases) / float(args.ws)


                        out1.write(scaff + "\t" +
                                   str(win_start) + "\t" +
                                   str(win_end) + "\t" +
                                   str(GD) + "\n")

                        num_CDS_bases = stop - win_end + 1
                        win_start = win_end
                        win_end += args.ws
                    elif start > win_end and scaff == oldscaff:
                        GD = float(num_CDS_bases) / float(args.ws)

                        out1.write(scaff + "\t" +
                                   str(win_start) + "\t" +
                                   str(win_end) + "\t" +
                                   str(GD) + "\n")

                        while start > win_end:
                            win_end += args.ws
                            win_start = win_end - args.ws
                            if start > win_end:
                                out1.write(scaff + "\t" +
                                           str(win_start) + "\t" +
                                           str(win_end) + "\t" +
                                           str(0.0) + "\n")

                        num_CDS_bases = 0

                        if all(k > win_start for k in [start, stop]) and all(k <= win_end for k in [start, stop]) and scaff == oldscaff:  # CDS is entirely within window bounds
                            num_CDS_bases += stop - start + 1
                        if start > win_start and start < win_end and stop > win_end and scaff == oldscaff:  # CDS is overhanging end of window
                            # Move onto next window, but parse overhang appropriately
                            num_CDS_bases += win_end - start + 1

                            GD = float(num_CDS_bases) / float(args.ws)

                            out1.write(scaff + "\t" +
                                       str(win_start) + "\t" +
                                       str(win_end) + "\t" +
                                       str(GD) + "\n")

                            num_CDS_bases = stop - win_end
                            win_start = win_end
                            win_end += args.ws

