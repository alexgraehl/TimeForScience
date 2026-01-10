// Density.java: By Charlie Vaske

// Topomap3D: A 3-D graph-viewer Cytoscape plug-in
// This source file is part of the Topomap3D plug-in.

// Copyright (C) 2005
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import java.awt.geom.Point2D;
import java.util.Random;


public class Density {
	public static final int MIN_PT = 0; // constant
	public static final int MAX_PT = 1; // constant
	
	private double[][] h; // height field

	private int totalSpaces;   // How many TOTAL grid spaces there are
	private int borderSpaces; // Out of TotalSpaces total spaces, how many spaces are border squares?

	public int getTotalSpaces()   { return totalSpaces; } // total spaces in EACH dimension
	public int getBorderSpaces() { return borderSpaces; } // border spaces in EACH dimension
	
	public double getNonBorderFraction() {
		// Returns the fraction of spaces in each linear dimension that are not border squares
		return ((double)totalSpaces - 2*borderSpaces) / (double)totalSpaces;
	}
	
	public Density(int xsize, int ysize) {
		h = new double[xsize][ysize];
	}
	
	public Density(Point2D[] points) {
		this(points, 3, 40, makeQuadraticKernel(3, -1, 0, 16));
	//this(points, 3, 250, makeQuadraticKernel(3, -1, 0, 16));
	}
	
	public Density(Point2D[] points, int theBorderSpaces, int theTotalSpaces, Density kernel) {
		this(theTotalSpaces, theTotalSpaces); // first: initialize a SQUARE array with some number of TotalSpaces
		
		this.totalSpaces  = theTotalSpaces; // TOTAL number of grid spaces
		this.borderSpaces = theBorderSpaces; // out of the total, how many are border spaces?
				
		// I think that borderSpaces is counted out of TotalSpaces from EACH direction
		// So there is a 2-unit border on each side of the gride (4 units total in each dimension) if you pass in borderSpaces = 2
		// So if TotalSpaces = 10, and borderSpaces = 2, then there are 6 non-border spaces in each dimension of the grid
		
		Point2D[] bounds = getBoundsOfPoints(points);
		
		if (points == null || kernel == null ||  bounds == null ||
			totalSpaces <= 0 || borderSpaces < 0 || totalSpaces <= 2*borderSpaces) { 
			System.err.println("Tried to create a density with an invalid parameter.");
			return;
		}
		
		double xrange = bounds[MAX_PT].getX() - bounds[MIN_PT].getX();
		double yrange = bounds[MAX_PT].getY() - bounds[MIN_PT].getY();
		
		double maxRange = Math.max(xrange, yrange);
		
		double spacesize = maxRange / (totalSpaces - 2*borderSpaces);
		
		double xbeg = bounds[MIN_PT].getX() - (borderSpaces*spacesize) - (maxRange-xrange)/2;
		double ybeg = bounds[MIN_PT].getY() - (borderSpaces*spacesize) - (maxRange-yrange)/2;
		
		int koffset = kernel.h.length / 2;
		
		for (int i = 0; i < points.length; ++i) {
			int xidx = (int)((points[i].getX() - xbeg) / spacesize);
			int yidx = (int)((points[i].getY() - ybeg) / spacesize);
			addDensity(kernel, xidx-koffset, yidx-koffset);
		}
	}
    
	static public Density createRandomDensity(int numRandomPoints) {
	  // Creates a random
		Random rands = new Random(); // default random seed is the current time in milliseconds
		
		Point2D[] randPoints = new Point2D[numRandomPoints];
		for (int i = 0; i < randPoints.length; i++) {
			randPoints[i] = new Point2D.Float(rands.nextFloat(), rands.nextFloat());
		}
		
		return new Density(randPoints);
	}
	
	public double[][] getHeightField() {
		return h;
	}
	
	public Density unsafeAddDensity(Density d, int x, int y) {
		for (int i = 0; i < d.h.length; ++i) {
			for (int j = 0; j < d.h[0].length; ++j) {
				h[x+i][y+j] += d.h[i][j];
			}
		}
		return this;
	}
	
	public Density addDensity(Density d, int x, int y) {
		int srci = 0;
		int srcj = 0;
		
		int dsti = x;
		int dstj = y;
		if (dsti < 0) {
			srci = -dsti;
			dsti = 0;
		}
		if (dstj < 0) {
			srcj = -dstj;
			dstj = 0;
		}
		
		int endi = d.h.length;
		int endj = d.h[0].length;
		if (dsti + endi > h.length) {
			endi = h.length - dsti;
		}
		if (dstj + endj > h[0].length) {
			endj = h[0].length - dstj;
		}
		
		for (; srci < endi; ++srci, ++dsti) {
			int oldsrcj = srcj;
			int olddstj = dstj;
			for (; srcj < endj; ++srcj, ++dstj) {
				h[dsti][dstj] += d.h[srci][srcj];
			}
			srcj = oldsrcj;
			dstj = olddstj;
		}
		
		return this;
	}
	
	public double maxValue() {
		// Returns the maximum value from the density
		double max = -Double.MAX_VALUE;
		for (int i = 0; i < h.length; i++) {
			for (int j = 0; j < h[0].length; j++) {
				if (h[i][j] > max) { max = h[i][j]; }
			}
		}
		return max;
	}
	
	public Density multiply(double d) {
		// Multiplies all elements by "d"
		if (h != null) {
			for (int i = 0; i < h.length; i++) {
				for (int j = 0; j < h[0].length; j++) {
					h[i][j] *= d;
				}
			}
		}
		return this;
	}
	
	public Density addScalar(double d) {
		// Adds "d" to all elements
		if (h != null) {
			for (int i = 0; i < h.length; i++) {
				for (int j = 0; j < h[i].length; j++) {
					h[i][j] += d;
				}
			}
		}
		return this;
	}
	
	static public Density makeQuadraticKernel(int radius, double a, double b, double c) {
		double[][] quadrant = new double[radius][radius];
		for (int i = 0; i < radius; ++i) {
			for (int j = 0; j < radius; ++j) {
				double r2 = i*i + j*j;
				double r = Math.sqrt(r2);
				quadrant[i][j] = a*r2 + b*r + c;
			}
		}
		int ksize = 2*radius - 1;
		Density r = new Density(ksize, ksize);
		
		for (int i = 0; i < radius; ++i) { 
			for (int j = 0; j < radius; ++j) {
				r.h[i][j] = quadrant[radius-i-1][radius-j-1];
			}
			for (int j = 1; j < radius; ++j) {
				r.h[i][radius+j-1] = quadrant[radius-i-1][j];
			}
		}
		for (int i = 1; i < radius; ++i) {
			for (int j = 0; j < radius; ++j) {
				r.h[radius+i-1][j] = quadrant[i][radius-j-1];
			}
			for (int j = 1; j < radius; ++j) {
				r.h[radius+i-1][radius+j-1] = quadrant[i][j];
			}
		}
		
		for (int i = 0; i < r.h.length; ++i) {
			for (int j = 0; j < r.h[i].length; ++j) {
				if (r.h[i][j] < 0.0) {r.h[i][j] = 0.0;}
			}
		}
		
		return (Density)r;
	} // end makeQuadraticKernel
	
	public static Point2D[] getBoundsOfPoints(Point2D[] nodes) {
		if (nodes == null || nodes.length == 0) { 
			System.err.println("Warning: getBoundsOfPoints(...) was called with a null argument.");
			return null;
		}
		
		double xmin = nodes[0].getX();
		double xmax = nodes[0].getX();
		double ymin = nodes[0].getY();
		double ymax = nodes[0].getY();
		
    // find bounds
		for (int i = 1; i < nodes.length; ++i) {
			double x = nodes[i].getX();
			double y = nodes[i].getY();
			
			if (x < xmin) {xmin = x;}
			if (x > xmax) {xmax = x;}
			if (y < ymin) {ymin = y;}
			if (y > ymax) {ymax = y;}
		}
		
		Point2D[] result = new Point2D[2];
		result[MIN_PT] = new Point2D.Double(xmin, ymin);
		result[MAX_PT] = new Point2D.Double(xmax, ymax);
		return result;
	} // end getBoundsOfPoints
	
	public String toString() {
		StringBuffer s = new StringBuffer();
		
		if (h == null) { return "Empty Density"; }
		
		for(int i = 0; i < h.length; ++i) {
			for (int j = 0; j < h[i].length; ++j) {
				s.append(h[i][j]).append("\t");
			}
			s.append("\n");
		}
		return s.toString();
	} // end toString
	
	public static void main(String[] argv) {
		Density qk = makeQuadraticKernel(3, -1, 0, 16);
		System.out.println("The kernel:\n" + qk);
		Density plains = new Density(25, 9);
		System.out.println("\n\nThe plain...\n" + plains);
		
		plains.addDensity(qk, 5, 2);
		System.out.println("\n\nAdding in the middle\n" + plains);
		
		System.out.println("\n\nAdding in upperleft corner\n" + 
						   (new Density(25,9)).addDensity(qk, 0, 0));
		
		System.out.println("\n\nAdding in upper right corner\n" + 
						   (new Density(25,9)).addDensity(qk, 0, 9-5));
		
		System.out.println("\n\nAdding in lower left corner\n" + 
						   (new Density(25,9)).addDensity(qk, 25-5, 0));
		
		System.out.println("\n\nAdding in lower right corner\n" + 
						   (new Density(25,9)).addDensity(qk, 25-5, 9-5));
		
		System.out.println("\n\nAdding overlapping upperleft corner\n" + 
						   (new Density(25,9)).addDensity(qk, -1, -1));
		
		System.out.println("\n\nAdding overlapping upper right corner\n" + 
						   (new Density(25,9)).addDensity(qk, -1, 9-5+1));
		
		System.out.println("\n\nAdding overlapping lower left corner\n" + 
						   (new Density(25,9)).addDensity(qk, 25-5+1, -1));
		
		System.out.println("\n\nAdding overlapping lower right corner\n" + 
						   (new Density(25,9)).addDensity(qk, 25-5+1, 9-5+1));
	} // end Density.main
	
} // class Density