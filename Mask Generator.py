import numpy as np
import gdspy #install gdspy
from gdspy import *
# https://gdspy.readthedocs.io/en/stable/geometry.html#classes
import matplotlib
matplotlib.use('TkAgg') #Prevent import error
import matplotlib.pyplot as plt
from functools import partial #This function just reduce the functionality of function
# https://docs.python.org/2/library/functools.html
import math

Cell = partial(Cell, exclude_from_current=True)
'''
Make a partial function from cell
Class in gdspy that you can look up at site-package
    exclude_from_current : bool
        If ``True``, the cell will not be automatically included in the
        current library.
'''

from text import Text as _Text
#Need to install freetype, no idea why we need this

import shapely
from shapely.geometry import Polygon as shapely_Polygon
from shapely.geometry import Point as shapely_Point
#Require installation of GEOS and then Shapely

import snakenmake
from snakenmake import *
#snake, profilometry_marks, alignment_cross, wafer, chip, outline

#Read the files of wires
Imposed_File = r'C:\Software\snakenbake-master\Transport_Cut.gds'
gdsii = gdspy.GdsLibrary()
file = gdsii.read_gds(Imposed_File)
TopLevelCell = file.top_level()[0]
Path = TopLevelCell.get_polygons() #From dxf to gds, it generate only all the lines, so I manually translate them into polygon with defined coordinate
for i in range(len(Path)): #Make them into a list, otherwise the index wont work
    Path[i] = Path[i].tolist()

List_Polygon, coordinatearray = [], None #Total polygon and temperary polygon
def Find_Next_Line(Path, LineEnd): #In the Path, find the line that can connect to LineEnd and return the line and its index in Path
    for line in Path:
        if LineEnd == line[0]:
            nextline = line[1]
            index = Path.index(line)
        elif LineEnd == line[1]:
            nextline = line[0]
            index = Path.index(line)
    return nextline, index
    
while Path != []: #Generate all the polygons
    if coordinatearray is None:
        coordinatearray = Path[0]
        Path.pop(0)
    else:
        if coordinatearray[-1] == coordinatearray[0]:
            List_Polygon.append(coordinatearray)
            coordinatearray = None
        else:
            newline, index = Find_Next_Line(Path, coordinatearray[-1])
            coordinatearray.append(newline)
            Path.pop(index)

for i in range(len(List_Polygon)): #Turn the polygon into real polygon
    List_Polygon[i] = shapely_Polygon(List_Polygon[i][0: -1])


#Parameters
FileDirectory = r'C:\Software\snakenbake-master\raymondtest.gds'
QRCode_Radius = 0.75 #Micron
QRCode_Spacing = 2 #Spacing between each Octagon representing the bit
QRCode_Separation = 30 #Micron, spacing between each QR code
QRCode_Sides = 8
DeviceCenter = (58550.0, 40039.0) #Micron
QRCode_CenterClearence = 150 #Micron
QRCode_CenterIndex = (49.5, 49.5)
QRCode_Range_X, QRCode_Range_Y = range(0, 100), range(0, 100)
Wire_Clearence = 40.0 #Micron


ALignmentMarker_Small_Spacing = 90.0
ALignmentMarker_Small_Width = 1.0
ALignmentMarker_Small_Length = 2.0

ALignmentMarker_Big_Spacing = 900.0
ALignmentMarker_Big_Width = 15.0
ALignmentMarker_Big_Length = 30.0

def plot_cell(cell, exclude=()):
    # FROM: https://github.com/heitzmann/gdspy/issues/42
    poly_dict = cell.get_polygons(by_spec=True)
    plt.figure(figsize=(20,12))
    for layer_datatype, polys in poly_dict.items():
        if layer_datatype[0] in exclude:
            continue
        for poly in polys:
            plt.fill(*poly.T, lw=0.5, ec='k', fc=(1,0,0,0.5))
    plt.axes().set_aspect('equal', 'datalim')
    
def write_gds_file(main_cell, filename, unit=1.0e-6, precision=1.0e-9):
    writer = gdspy.GdsWriter(filename, unit=unit, precision=precision)
    for cell in [main_cell] + list(main_cell.get_dependencies(True)):
        writer.write_cell(cell)
    writer.close() 
    
def Octagon(center, radius, layer = 10): #This generate a octagon that will be used as a single qr code
    points = [(center[0] + radius * np.sin(np.pi / 8), center[1] - radius * np.cos(np.pi / 8)), 
              (center[0] + radius * np.cos(np.pi / 8), center[1] - radius * np.sin(np.pi / 8)), 
              (center[0] + radius * np.cos(np.pi / 8), center[1] + radius * np.sin(np.pi / 8)), 
              (center[0] + radius * np.sin(np.pi / 8), center[1] + radius * np.cos(np.pi / 8)), 
              (center[0] - radius * np.sin(np.pi / 8), center[1] + radius * np.cos(np.pi / 8)), 
              (center[0] - radius * np.cos(np.pi / 8), center[1] + radius * np.sin(np.pi / 8)), 
              (center[0] - radius * np.cos(np.pi / 8), center[1] - radius * np.sin(np.pi / 8)), 
              (center[0] - radius * np.sin(np.pi / 8), center[1] - radius * np.cos(np.pi / 8))]
    return Polygon(points, layer)
    
def Multigon(center, radius, sides, layer = 10):
    points = []
    for i in range(sides):
        points.append((center[0] + radius * np.sin(((2 * i) + 1) * np.pi / sides), center[1] + radius * np.cos(((2 * i) + 1) * np.pi / sides)))
    return Polygon(points, layer)

def Add_QRCode(Cell, center, radius, spacing, xIndex, yIndex, layer = 10):
    xCoordinate, yCoordinate = center[0], center[1]
    
    #Generate alignment QR codes
    Cell.add(Multigon((xCoordinate - 2.5 * spacing, yCoordinate - 2.5 * spacing), radius, QRCode_Sides))
    Cell.add(Multigon((xCoordinate - 1.5 * spacing, yCoordinate - 2.5 * spacing), radius, QRCode_Sides))
    Cell.add(Multigon((xCoordinate + 2.5 * spacing, yCoordinate - 2.5 * spacing), radius, QRCode_Sides))

    #Generate x portion of our QR codes
    binary_x = GenerateBinaray(xIndex)
    for i in range(8):
        if binary_x[i] == '1':
            dot_position_x, dot_position_y = xCoordinate - 1.5 * spacing + spacing * (i % 4), yCoordinate + 1.5 * spacing - spacing * int(i / 4)
            Cell.add(Multigon((dot_position_x, dot_position_y), radius, QRCode_Sides))
            
    #Generate y portion of our QR codes
    binary_y = GenerateBinaray(yIndex)
    for i in range(8):
        if binary_y[i] == '1':
            dot_position_x, dot_position_y = xCoordinate - 1.5 * spacing + spacing * (i % 4), yCoordinate - 0.5 * spacing - spacing * int(i / 4)
            Cell.add(Multigon((dot_position_x, dot_position_y), radius, QRCode_Sides))
    
def Add_AlginmentMarker(Cell, center, spacing, width, length, layer = 10):
    Add_Cross(Cell, (center[0] + spacing / 2, center[0] + spacing / 2), width, length, layer)
    Add_Cross(Cell, (center[0] - spacing / 2, center[0] + spacing / 2), width, length, layer)
    Add_Cross(Cell, (center[0] + spacing / 2, center[0] - spacing / 2), width, length, layer)
    Add_Cross(Cell, (center[0] - spacing / 2, center[0] - spacing / 2), width, length, layer)
    
def Add_Cross(Cell, center, width, length, layer = 10):
    Cell.add(Rectangle((center[0], center[1]),(center[0] + width, center[1] + length),layer=10))
    Cell.add(Rectangle((center[0] + width, center[1]),(center[0] + length, center[1] + width),layer=10))
    Cell.add(Rectangle((center[0], center[1]),(center[0] - width, center[1] - length),layer=10))
    Cell.add(Rectangle((center[0] - width, center[1]),(center[0] - length, center[1] - width),layer=10))
    
def GenerateBinaray(integervalue):
    binary = bin(integervalue)[2:] #Generate binary code as string
    while len(binary) < 8: #Supplement 0 to fit 8 bit
        binary = '0' + binary
    return binary

def plot_cell(cell, exclude=()):
    # FROM: https://github.com/heitzmann/gdspy/issues/42
    poly_dict = cell.get_polygons(by_spec=True)
    plt.figure(figsize=(20,12))
    for layer_datatype, polys in poly_dict.items():
        if layer_datatype[0] in exclude:
            continue
        for poly in polys:
            plt.fill(*poly.T, lw=0.5, ec='k', fc=(1,0,0,0.5))
    plt.axes().set_aspect('equal', 'datalim')
    
def write_gds_file(main_cell, filename, unit=1.0e-6, precision=1.0e-9):
    writer = gdspy.GdsWriter(filename, unit=unit, precision=precision)
    for cell in [main_cell] + list(main_cell.get_dependencies(True)):
        writer.write_cell(cell)
    writer.close()

#Actual Code of drawing start here

main_cell = Cell('main')

new_cell = Cell('new')

for x in QRCode_Range_X:
    for y in QRCode_Range_Y:
        QRCode_PositionX, QRCode_PositionY = DeviceCenter[0] + (x - QRCode_CenterIndex[0]) * QRCode_Separation, DeviceCenter[1] + (y - QRCode_CenterIndex[1]) * QRCode_Separation
        add_flag = True
        if math.sqrt((QRCode_PositionX - DeviceCenter[0]) ** 2 + (QRCode_PositionY - DeviceCenter[1]) ** 2) <= QRCode_CenterClearence:
            add_flag = False
        if add_flag:
            point = shapely_Point(QRCode_PositionX, QRCode_PositionY)
            for i in range(len(List_Polygon)):
                if add_flag and List_Polygon[i].distance(point) <= Wire_Clearence:
                    add_flag = False
        if add_flag:
            Add_QRCode(new_cell, (QRCode_PositionX, QRCode_PositionY), QRCode_Radius, QRCode_Spacing, x, y, layer = 10)

# Add_AlginmentMarker(new_cell, DeviceCenter, ALignmentMarker_Small_Spacing, ALignmentMarker_Small_Width, ALignmentMarker_Small_Length, layer = 10)

# Add_AlginmentMarker(new_cell, DeviceCenter, ALignmentMarker_Big_Spacing, ALignmentMarker_Big_Width, ALignmentMarker_Big_Length, layer = 10)

#There is no code in this cell, this is just a reminder that adding the nanosquid logo would 
#be cool, and you should do that in L-Edit

main_cell.add(CellReference(new_cell,(0,0)))

write_gds_file(main_cell, FileDirectory)
