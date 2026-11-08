#! /usr/bin/env python3
# Author: Izaak Neutelings (August 2025)
# Description:
#   Create a Venn diagram as a TH2Poly that can be filled.
import math
from math import pi, sqrt, cos, sin, atan2, acos, degrees, ceil
from array import array
from enum import IntEnum
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True # to avoid conflict with argparse
from ROOT import gROOT, gStyle, TCanvas, TH2Poly, TGraph, TEllipse, TLatex,\
                 kBlack, kRed, kBlue
gROOT.SetBatch(True) # do not open GUI windows


VennBins = IntEnum('VennBins', [
  ('A', 1),('B', 2),('C', 3),
  ('AnB', 4),('AnC', 5),('BnC', 6),
  ('AnBnC', 7)
])
vbincenters_dict = { }


def TH2Venn(name='venn', title='', #'venn;x;y;',
            xmin=None, xmax=None, ymin=None, ymax=None,
            r0=None, R=1, nsteps=200, marg=0.25, rotate=0):
  """Create Venn diagram of three intersecting circles as TH2Poly object."""
  if r0 is None:
    r0 = 0.7*R # radial coordinate of centers
  O1 = polar2xy(r0,rotate+150)
  O2 = polar2xy(r0,rotate+30)
  O3 = polar2xy(r0,rotate-90)
  centers = [O1,O2,O3]
  if xmin is None:
    xmin = min(o[0] for o in centers)-(1+marg)*R
  if xmax is None:
    xmax = max(o[0] for o in centers)+(1+marg)*R
  if ymin is None:
    ymin = min(o[1] for o in centers)-(1+marg)*R
  if ymax is None:
    ymax = max(o[1] for o in centers)+(1+marg)*R
  hist = TH2Poly(name,title,xmin,xmax,ymin,ymax)
  hist.SetFloat(True)  # allow non-integer bin contents
  hist.SetStats(False) # do not draw statistics box
  
  # FULL CIRCLE
  dt = 2*pi/nsteps # steps width of angle
  for i, o1 in enumerate(centers):
    o2 = centers[(i+1)%len(centers)]
    o3 = centers[(i+2)%len(centers)]
    x, y = [ ], [ ]
    _, t1 = getintersectionangles(*o1,*o2,R)
    t2, _ = getintersectionangles(*o1,*o3,R)
    _, t3 = getintersectionangles(*o3,*o1,R)
    _, t4 = getintersectionangles(*o3,*o2,R)
    t5, _ = getintersectionangles(*o2,*o3,R)
    t6, _ = getintersectionangles(*o2,*o1,R)
    getxyarc(*o1,R,t1=t1,t2=t2,dt=dt,x=x,y=y)
    getxyarc(*o3,R,t1=t3,t2=t4,dt=dt,x=x,y=y,s=-1)
    getxyarc(*o2,R,t1=t5,t2=t6,dt=dt,x=x,y=y,s=-1)
    graph = TGraph(len(x), array('d',x), array('d',y))
    graph.SetLineWidth(1)
    graph.SetLineColor(kBlack)
    hist.AddBin(graph)
  
  # INTERSECTION between two circles
  centers = [O1,O3,O2]
  for i, o1 in enumerate(centers):
    o2 = centers[(i+2)%len(centers)]
    o3 = centers[(i+1)%len(centers)]
    x, y = [ ], [ ]
    _, t1 = getintersectionangles(*o1,*o2,R)
    _, t2 = getintersectionangles(*o1,*o3,R)
    t3, _ = getintersectionangles(*o3,*o1,R)
    _, t4 = getintersectionangles(*o3,*o2,R)
    t5, _ = getintersectionangles(*o2,*o3,R)
    t6, _ = getintersectionangles(*o2,*o1,R)
    getxyarc(*o1,R,t1=t1,t2=t2,dt=dt,x=x,y=y,s=-1)
    getxyarc(*o3,R,t1=t3,t2=t4,dt=dt,x=x,y=y)
    getxyarc(*o2,R,t1=t5,t2=t6,dt=dt,x=x,y=y,s=-1)
    graph = TGraph(len(x), array('d',x), array('d',y))
    graph.SetLineWidth(1)
    graph.SetLineColor(kBlack)
    hist.AddBin(graph)
  
  # INTERSECTION between three circles
  x, y = [ ], [ ]
  t1, _ = getintersectionangles(*O1,*O2,R)
  _, t2 = getintersectionangles(*O1,*O3,R)
  t3, _ = getintersectionangles(*O3,*O1,R)
  _, t4 = getintersectionangles(*O3,*O2,R)
  t5, _ = getintersectionangles(*O2,*O3,R)
  _, t6 = getintersectionangles(*O2,*O1,R)
  getxyarc(*O1,R,t1=t1,t2=t2,dt=dt,x=x,y=y)
  getxyarc(*O3,R,t1=t3,t2=t4,dt=dt,x=x,y=y)
  getxyarc(*O2,R,t1=t5,t2=t6,dt=dt,x=x,y=y)
  graph = TGraph(len(x), array('d',x), array('d',y))
  graph.SetLineWidth(1)
  graph.SetLineColor(kBlack)
  hist.AddBin(graph)
  
  hist.SetMarkerSize(2)
  return hist
  

def getxyarc(Ox,Oy,R,t1=0,t2=2*pi,dt=0.025,s=1,x=None,y=None):
  """
  Compute list of (x,y) coordinates following an arc
  of circle with center (Ox,Oy) and radius R
  between t1 and t2 (radians) on a circle.
  Default counterclockwise (s=1), otherwise clockwise (s=-1).
  """
  if x is None:
    x = [ ]
  if y is None:
    y = [ ]
  t = t1 # start
  while s*t2<s*t1: # fold
    t2 += s*2*pi
  dt = (t2-t1)/round(s*(t2-t1)/dt) # shrink/stretch
  while abs((t-t2)%(2*pi))>=abs(dt):
    x.append(Ox+R*cos(t))
    y.append(Oy+R*sin(t))
    t += dt
  return x, y
  

def getvbincenters(R=1, r0=None, rotate=0):
  global vbincenters_dict
  if r0 is None:
    r0 = 0.7*R # radial coordinate of centers
  vbincenters_dict = {
    VennBins.AnBnC:  polar2xy(0.05*r0,      rotate+90),
    VennBins.AnB:    polar2xy(0.9*r0+0.02*R,rotate+90),
    VennBins.AnC:    polar2xy(0.9*r0+0.02*R,rotate-150),
    VennBins.BnC:    polar2xy(0.9*r0+0.02*R,rotate-30),
    VennBins.A:      polar2xy(    r0+0.35*R,rotate+150),
    VennBins.B:      polar2xy(    r0+0.35*R,rotate+30),
    VennBins.C:      polar2xy(    r0+0.35*R,rotate-90),
  }
  return vbincenters_dict
  

#def getangle(x1,y1,x2,y2):
#  return atan2(y2-y1,x2-x1) # radians
  

#def angle_diff(a,b):
#  return (a-b)%(2*pi)
  

def polar2xy(r,theta,r0=(0,0)):
  """Convert polar (r,theta) coordinates to (x,y) coordinates"""
  theta = pi/180*theta
  return r0[0]+r*math.cos(theta), r0[1]+r*math.sin(theta)
  

def getintersectionangles(x1,y1,x2,y2,R):
  """Compute intersection points of two circles with radii R,
  and centers O1=(x1,y1), O2=(x1,y1). Return the angles between
  the horizontal and O1 to the intersection points."""
  dx, dy = (x2-x1), (y2-y1)
  dang = acos(sqrt(dx**2+dy**2)/2/R)
  ang  = atan2(dy, dx) # angle (O1)-(O2) w.r.t. horizontal
  ang1 = ang - dang
  ang2 = ang + dang
  #print(f">>> getintersectionangles: O1=({x1:.2f},{y1:.2f}), O2=({x2:.2f},{y2:.2f}) "
  #      f" => ang={degrees(ang)}, dang={degrees(dang)}")
  return ang1, ang2
  

def getintersection(x1,y1,x2,y2,r):
  """Get intersection points of two overlapping circles."""
  dx, dy = (x2-x1), (y2-y1)
  d  = sqrt(dx**2 + dy**2)
  h  = sqrt((r/d)**2-1/4)
  pin  = ( (x1+x2)/2+h*dy, (y1+y2)/2-h*dx )
  pout = ( (x1+x2)/2-h*dy, (y1+y2)/2+h*dx )
  return pin, pout
  

def drawvenncircles(R=1,r0=None,lwidth=7,colors=[kBlack,kBlue,kRed]):
  """Draw Venn diagram circles"""
  if r0 is None:
    r0 = 0.7*R # radial coordinate of centers
  O1 = polar2xy(r0,150)
  O2 = polar2xy(r0, 30)
  O3 = polar2xy(r0,-90)
  circles = [ ]
  for i, o in enumerate([O1,O2,O3]):
    circle = TEllipse(*o,R)
    circle.SetFillColorAlpha(0,0)
    circle.SetLineWidth(lwidth)
    if colors:
      circle.SetLineColor(colors[i])
    circle.Draw('SAME L')
    circles.append(circle)
  return circles # avoid garbage collection
  

def drawpoint(x1,y1,x2,y2,r,n=1):
  pin, pout = getintersection(x1,y1,x2,y2,r)
  graphs = [ ]
  texts = [ ]
  for i, p in [(-1,pin),(1,pout)]:
    graph = TEllipse(*p,0.018)
    graph.SetLineWidth(2)
    graph.SetLineColor(kBlack)
    graph.SetFillColor(kBlack)
    graph.Draw('SAME F')
    graphs.append(graph)
    latex = TLatex()
    latex.SetTextSize(0.025)
    latex.SetTextAlign(22)
    latex.SetTextFont(42)
    latex.SetTextColor(kRed)
    latex.DrawLatex(*p,f"{i*n}")
    texts.append(latex)
  return graphs, texts
  

def getvenncanvas(cname='canvas', width=600, height=570, drawpoints=False, hideaxes=True,
                  format=".0f", R=1, **kwargs):
  """Return canvas and Venn histogram."""
  
  # CANVAS
  #xmin, xmax = -0.5, 0.5
  canvas = TCanvas(cname, cname, width, height)
  #frame = canvas.DrawFrame(xmin, xmin, xmax, xmax)
  #frame.GetXaxis().SetRangeUser(-1,1)
  #frame.GetYaxis().SetRangeUser(-1,1)
  #frame.Draw('AXIS')
  
  # VENN
  if hideaxes:
    kwargs.setdefault('marg',0.3)
  hist = TH2Venn(R=R,**kwargs)
  hist.SetStats(0)
  #hist.SetMinimum(0)
  #hist.SetMaximum(hist.GetNumberOfBins())
  
  # DRAW
  gStyle.SetPaintTextFormat(format)
  hist.SetMarkerSize(2)
  hist.SetMarkerColor(kRed)
  if hideaxes:
    canvas.SetMargin(0,0,0,0) # LRBT
    canvas.SetFrameLineColor(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetTicks(False,False)
    for axis in [hist.GetXaxis(),hist.GetYaxis()]:
      axis.SetLabelSize(0)
      axis.SetTickLength(0)
      axis.SetAxisColor(0,0) # transparent
  else:
    canvas.SetMargin(0.11,0.08,0.10,0.08) # LRBT
  
  # DRAW INTERSECTION POINTS (for debugging)
  objects = [ ] # to prevent garbage collection
  if drawpoints:
    r0 = kwargs.get('r0',0.7*R) # radial coordinate of center
    O1 = polar2xy(r0, 30)
    O2 = polar2xy(r0,150)
    O3 = polar2xy(r0,-90)
    c1 = drawpoint(*O2,*O1,R,1)
    c3 = drawpoint(*O3,*O2,R,2)
    c2 = drawpoint(*O1,*O3,R,3)
  
  canvas.Update()
  if drawpoints:
    return canvas, hist, objects
  return canvas, hist
  

if __name__ == '__main__':
  gStyle.SetPaintTextFormat(".0f")
  
  # WITHOUT AXES
  canvas, hist, objects = getvenncanvas(drawpoints=True,hideaxes=False)
  for i in range(1,hist.GetNumberOfBins()+1):
    hist.SetBinContent(i,i)
  hist.Draw('COLZ TEXT L') 
  canvas.SaveAs("vennhist_axes.png")
  canvas.SaveAs("vennhist_axes.pdf")
  canvas.Close()
  
  # WITH AXES
  canvas, hist, objects = getvenncanvas(drawpoints=True,hideaxes=True)
  for i in range(1,hist.GetNumberOfBins()+1):
    hist.SetBinContent(i,i)
  hist.Draw('COL TEXT L')
  canvas.SaveAs("vennhist.png")
  canvas.SaveAs("vennhist.pdf")
  canvas.Close()
  
