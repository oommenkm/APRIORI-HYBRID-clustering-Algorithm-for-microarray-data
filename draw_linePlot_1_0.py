#!/usr/bin/env python
from reportlab.graphics.shapes import Drawing
from reportlab.lib import colors
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.graphics.charts.legends import Legend
from reportlab.graphics.charts.textlabels import Label

d = Drawing(620, 450)

lp = LinePlot()
lp.width = 520
lp.height = 320
lp.x = 70
lp.y = 100

FH = open("finalOutPut.txt")

data_dict = eval(FH.read())

gcs = []
for g in data_dict.keys():
    gcs.append(int(g))
    
gcs.sort()
gc_strs = []

pos_pers = []
neg_pers = []
nor_pos_pers = []
nor_neg_pers = []

for gc in gcs:
    gn = str(gc)
    gc_strs.append(gn)
    pos_pers.append((gc, round(data_dict[gn]['pos_per'], 2)))
    neg_pers.append((gc, round(data_dict[gn]['neg_per'], 2)))
    nor_pos_pers.append((gc, round(35.483870967741936, 2)))
    nor_neg_pers.append((gc, round(64.516129032258064, 2)))
    print gn, round(data_dict[gn]['pos_per']), round(data_dict[gn]['neg_per'])

lp.data = [pos_pers, neg_pers, nor_pos_pers, nor_neg_pers]
print lp.data

lp.joinedLines = 1
lp.lines.symbol = makeMarker('Circle')
lp.lineLabelFormat = '%2.2f'
lp.strokeColor = colors.black
lp.xValueAxis.valueMin = 0
lp.xValueAxis.valueMax = 2100
lp.xValueAxis.labelTextFormat = '%2.0f'
lp.yValueAxis.valueMin = 0
lp.yValueAxis.valueMax = 100
lp.yValueAxis.valueStep = 10

xlbl = Label()
xlbl.setText("No. of Genes")
xlbl.setOrigin(310, 72)

ylbl = Label()
ylbl.setText("Percentage\n      (%)")
ylbl.setOrigin(28, 260)

lp.lines[0].strokeColor = colors.darkgreen
lp.lineLabels[0].strokeColor = colors.darkgreen
lp.lines[1].strokeColor = colors.tomato
lp.lineLabels[1].strokeColor = colors.tomato
lp.lines[2].strokeColor = colors.aquamarine
lp.lineLabels[2].strokeColor = colors.aquamarine
lp.lines[3].strokeColor = colors.purple
lp.lineLabels[3].strokeColor = colors.purple

lgnd = Legend()
lgnd.x = 70
lgnd.y = 63
lgnd.autoXPadding = 45
lgnd.colorNamePairs = [ (lp.lines[0].strokeColor, 'Positive'),
                        (lp.lines[1].strokeColor, 'Negative'),
                        (lp.lines[2].strokeColor, 'Norm. Pos.'),
                        (lp.lines[3].strokeColor, 'Norm. Neg.')
                      ]
lgnd.demo()

d.add(lp)
d.add(lgnd)
d.add(xlbl)
d.add(ylbl)
d.save(fnRoot='testLinePlot1', formats=['png', 'pdf'])