add_library('pdf')

# The inputs for this script are lists of the outputs from the CSI.py script
# run using all species as the reference.
hii_list = loadStrings("./hii.txt")
cicada_list = loadStrings("./cicada_hii.txt")

out = createWriter("8_22_17.txt")
beginRecord(PDF, "8_22_17.pdf")
size(1000, 500)
background(255)

# These set the colors for the heatmap, corresponding to the CSI score
color_0 = 0xFF087FDB
color_50 = 0xFFFFF400
color_100 = 0xFFFF0000

txt_dict = {}
sp_list = []
# This makes a list of all species, and the cicada CSI files corresponding to
# them.
for x in range(0, len(cicada_list)):
    temp = cicada_list[x].split("\t")
    txt_dict[temp[0]] = temp[1]
    sp_list.append(temp[0])

cicada_dict = {}
# This opens the cicada CSI files and reads the values into a dictionary.
for ref in sp_list:
    cicada_dict[ref] = {}
    file = loadStrings(txt_dict[ref])
    for x in range(1, len(file)):
        temp = file[x].split("\t")
        target = temp[1]
        hii = float(temp[2])
        cicada_dict[ref][target] = hii
    cicada_dict[ref][ref] = float(1)

# Defines the size of individual cells of the heatmap, and where they'll start.
cell_width = 50
cell_height = 50
x_start = 10
y_start = 10
y = 0

# Makes the actual heatmap for the cicada CSI. Also writes a file with the
# corresponding values.
for sp in sp_list:
    x = 0
    out.print("%s\t" % sp)
    for target in sp_list:
        if cicada_dict[sp][target] <= 0.5:
            fill_color = lerpColor(color_0, color_50, cicada_dict[sp][target] / 0.5)
            fill(fill_color)
            stroke(fill_color)
        else:
            fill_color = lerpColor(color_50, color_100, (cicada_dict[sp][target] - 0.5) / 0.5)
            fill(fill_color)
            stroke(fill_color)
        rect(x_start + x, y_start + y, cell_width, cell_height)
        out.print("%s\t" % cicada_dict[sp][target])
        x += cell_width
    y += cell_height
    out.print("\n")

txt_dict = {}
# Makes a dictionary with the CSI files for each species.
for x in range(0, len(hii_list)):
    temp = hii_list[x].split("\t")
    txt_dict[temp[0]] = temp[1]

hii_dict = {}
# Opens the CSI files and inputs their values to a dictionary.
for ref in sp_list:
    hii_dict[ref] = {}
    file = loadStrings(txt_dict[ref])
    for x in range(1, len(file)):
        temp = file[x].split("\t")
        target = temp[2]
        hii = float(temp[4])
        if target not in hii_dict[ref]:
            hii_dict[ref][target] = []
        if ref == target:
            hii_dict[ref][target].append(float(1))
        else:
            hii_dict[ref][target].append(hii)

x_start += x_start + (cell_width * len(sp_list))
print x_start, len(sp_list)
y = 0

out.print("\n\n")

# Makes the heatmap for the regular CSI values.
for sp in sp_list:
    x = 0
    out.print("%s\t" % sp)
    for target in sp_list:
        avg = sum(hii_dict[sp][target]) / len(hii_dict[sp][target])
        print sp, target, avg
        if avg <= 0.5:
            fill_color = lerpColor(color_0, color_50, avg / 0.5)
            fill(fill_color)
            stroke(fill_color)
        else:
            fill_color = lerpColor(color_50, color_100, (avg - 0.5) / 0.5)
            fill(fill_color)
            stroke(fill_color)
        rect(x_start + x, y_start + y, cell_width, cell_height)
        out.print("%s\t" % avg)
        x += cell_width
    y += cell_height
    out.print("\n")

x_start += 10 + (cell_width * len(sp_list))
print x_start, len(sp_list)

strokeWeight(1.5)
# The following two for loops make the scale bar on the right of the figure.
for y in range(0, (cell_height * len(sp_list) / 2)):
    inter = float(y) / ((cell_height * len(sp_list)) / 2)
    stroke_color = lerpColor(color_0, color_50, inter)
    stroke(stroke_color)
    line(x_start + 20, y_start + (cell_height * len(sp_list)) - y, x_start + 40, y_start + (cell_height * len(sp_list)) - y)

for y in range(0, (cell_height * len(sp_list) / 2)):
    inter = float(y) / ((cell_height * len(sp_list)) / 2)
    stroke_color = lerpColor(color_50, color_100, inter)
    stroke(stroke_color)
    line(x_start + 20, y_start + ((cell_height * len(sp_list)) / 2) - y, x_start + 40, y_start + ((cell_height * len(sp_list)) / 2) - y)


endRecord()
out.flush()
out.close()
print "All done"
