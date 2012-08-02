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

//
//  LandscapeShape3D.java
//  Topomap3D
//
//  Created by alexgraehl on 9/22/05.

// windowing stuff
import java.awt.event.*;
import java.awt.GraphicsConfiguration;
import javax.swing.JOptionPane;

// Java stuff
import javax.vecmath.*;
import java.lang.Math;
import java.util.*;
import java.awt.geom.Point2D;

// 3d things
import com.sun.j3d.utils.applet.MainFrame; 
import com.sun.j3d.utils.universe.*;
import javax.media.j3d.*;
import com.sun.j3d.utils.geometry.*; // Sphere, Cone, Text2D, etc.

// GINY graph library
import giny.model.*;
import giny.view.GraphView;
import giny.view.EdgeView;
import giny.view.NodeView;

// Cytoscape stuff
import cytoscape.util.*;
import cytoscape.view.*;
import cytoscape.plugin.CytoscapePlugin;
import cytoscape.util.CytoscapeAction;
import cytoscape.Cytoscape;
import cytoscape.CyNetwork;
import cytoscape.CyNode;
import cytoscape.CyEdge;
import cytoscape.view.CyNetworkView;
import cytoscape.data.Semantics;

// http://download.java.net/media/java3d/javadoc/1.3.2/javax/media/j3d/TriangleStripArray.html
// http://download.java.net/media/java3d/javadoc/1.4.0-latest/com/sun/j3d/utils/behaviors/mouse/MouseBehavior.html

public class LandscapeShape3D extends Shape3D {
	
	// The LandscapeShape3D shape is defined in the unit cube from (0, 0, 0) to (1, 1, 1).
	// It is normally DRAWN in the cube from (-0.5, -0.5, 0) to (+0.5, +0.5, 1) (so the X-Y intersection is in the dead center of the view, and the terrain elevates up the Z axis)
	
	int gridVerticesX; // instance variables
	int gridVerticesY;
	
	float maxTerrainValue;
	float[][] terrain;   // local variable
	Density   dens;      // we need to remember the density that was passed in
	
	LandscapeShape3D(Density d) { // Constructor
		this.dens = d;
		this.maxTerrainValue = Float.MIN_VALUE; // maximum terrain value, later used to calculate height and colors
		this.setTerrainFromDensity(d);
		this.setGeometry(makeLandscapeGeometry());
		this.setAppearance(makeLandscapeAppearance());
		
		this.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE); // Note that only one capability bit may be set per 
		//this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);  // method invocation--capability bits cannot be ORed together.
		this.setCapability(Shape3D.ALLOW_APPEARANCE_WRITE);  // method invocation--capability bits cannot be ORed together.
		
	} // end of LandscapeShape3D constructor 
	
	public void randomizeTerrain() {
		this.setTerrainFromDensity(null);
	}
	
	public void setTerrainFromDensity(Density d) {
		// Changed non-local variables: terrain, gridVerticesX, gridVerticesY, maxTerrainValue
		
		if (d == null || d.getHeightField() == null) { // fall back on the random terrain if there is no specified density
			System.err.println("Alex: WARNING: no specified density file, or else empty density. Using randomized density instead.");
			d = Density.createRandomDensity(100);
		}
		
		double h[][] = d.getHeightField();
		
		assert (h != null) : "h was null in setTerrainFromDensity";
		
		this.gridVerticesX = h.length;
		this.gridVerticesY = h[0].length;
		
		this.terrain = new float[gridVerticesX][gridVerticesY];
		for(int i = 0; i < terrain.length; ++i) {
			for (int j = 0; j < terrain[i].length; ++j) {
				// We flip the Y coordinate here, because Cytoscape uses a different coordinate system
				// from what we like in the graphics world. In Cytoscape, the top left is (0, 0), and Y increases
				// going down the screen. In our 3D system, Y increases going up the screen, like a normal graph.
				// For display, we flip it so that we match Cytoscape's coordinate system.
				//  /\                                                    +---> x
				// y|        <--This is what we like.                    y|
				//  +----> x             This is how Cytoscape works.-->  \/
				int flippedJ = (this.gridVerticesY - j - 1);
				terrain[i][j] = (float)h[i][flippedJ];
				if (terrain[i][j] > this.maxTerrainValue) { this.maxTerrainValue = terrain[i][j];}
			}
		}
	}
	
	public Geometry makeLandscapeGeometry() {
		if (null == this.terrain) {
			System.out.println("makeLandscapeShape3D: attempting to operate on a not-yet-initialized terrain array. Initializing with a random terrain.");
			this.setTerrainFromDensity(null);
		}
		
		// Each strip is a vertical strip (i.e., parallel to the Y-axis) of triangles.
		// Strips are placed next to each other along the X axis
		int vertexCountPerStrip = (gridVerticesY) * 2;
		int totalVertexCounts   = vertexCountPerStrip * (gridVerticesX - 1);
		
		//System.out.println("Testing..." + vertexCountPerStrip);
		
		TriangleStripArray tris = (TriangleStripArray) this.getGeometry();
		if (null == tris) {
			// If the triangle strip array is not ALREADY set, then make a new one
			int triangleStripCounts[] = new int[ this.gridVerticesX - 1 ];
			for (int i = 0; i < (this.gridVerticesX - 1); i++) {
				triangleStripCounts[i] = vertexCountPerStrip;
			}
			
			tris = new TriangleStripArray(totalVertexCounts,
										  TriangleArray.COORDINATES | TriangleArray.COLOR_3,
										  triangleStripCounts);
			tris.setCapability(GeometryArray.ALLOW_COORDINATE_WRITE);
			tris.setCapability(GeometryArray.ALLOW_COORDINATE_READ);
			tris.setCapability(GeometryArray.ALLOW_COLOR_READ);  // Only one capability bit may be set per 
			tris.setCapability(GeometryArray.ALLOW_COLOR_WRITE); // invocation--bits cannot be ORed together.
			//tris.setCapability(GeometryArray.ALLOW_NORMAL_READ);
			//tris.setCapability(GeometryArray.ALLOW_NORMAL_WRITE);
			tris.setCapability(GeometryArray.ALLOW_COUNT_READ);
		}
		
		//System.out.println("Vertex count is: " + tris.getValidVertexCount() + " and also " + tris.getVertexCount());
		
		if (gridVerticesX < 2 || gridVerticesY < 2) {
			System.err.println("gridVerticesX or gridVerticesY was too small. X: " + gridVerticesX + ", Y: " + gridVerticesY);
			System.exit(1);
		}
		
		float yStep = 1.0f / (float)this.gridVerticesY; // total size is totalSize (no matter what)
		float xStep = 1.0f / (float)this.gridVerticesX; // total size is totalSize (no matter what)
		float zScale = 1.0f / this.maxTerrainValue; // sort of ad-hoc, but scales so that the maximum peak is always the same height
		
		//int counter = 0;
		float epsilon = 0.05f * this.maxTerrainValue;
		float epsmax = this.maxTerrainValue + epsilon;
		
		for (int indexInTriArray=0, x=0; x < (this.gridVerticesX-1); x++) { // note the "-1" in the condition here, but NOT in the y loop below!
			for (int y = 0; y < this.gridVerticesY; y++) {
				// y runs one more time than you might expect, because we have to finish off the top of each triangle
				// (points 2 and 3 in the diagram below)
				//   2-------3  For a 1x1 grid, we draw two triangles, in the
				//   | \     |  order shown left.
				//   |   \   |  Keep extending the grid vertically up in order to add more grid units.
				//   |     \ |
				//   0-------1
				
				// Set the LEFT point (0 on the diagram above)
				float intensity1 = (this.terrain[x][y] + epsilon) / epsmax; // ranges from epsilon to 1
				tris.setCoordinate(indexInTriArray, new Point3f(x*xStep,
																y*yStep, 
																zScale * this.terrain[x][y]));
				tris.setColor(indexInTriArray, colorFromIntensity(intensity1));
				
				// Set the RIGHT point (1 on the diagram above)
				float intensity2 = (this.terrain[x+1][y] + epsilon) / epsmax; // ranges from epsilon to 1
				tris.setCoordinate(indexInTriArray+1, new Point3f((x+1)*xStep,
																  y*yStep,
																  zScale * this.terrain[x+1][y]));
				tris.setColor(indexInTriArray+1, colorFromIntensity(intensity2));
				
				// y   P
				// |
				// |   R    Y    // Testing the orientation
				// + -=-=-=-=-
				//      x
				//if (x == 5 && y == 5)  tris.setColor(indexInTriArray, new Color3f(1.0f, 0, 0)); 
				//if (x == 10 && y == 5) tris.setColor(indexInTriArray, new Color3f(1.0f, 1.0f, 0));
				//if (x == 5 && y == 10) tris.setColor(indexInTriArray, new Color3f(1.0f, 0, 1.0f));
				
				indexInTriArray += 2;
				
				//System.err.println( "This is " + (x+1)*xStep + ", " + y*yStep + ".");
			}
		}
		return tris;
	}
	
	public static Color3f colorFromIntensity(float intensity) {
		// intensity must range from 0.0f to 1.0f (inclusive: [0..1])
		return new Color3f(
						   (float) // RED
						   Math.max(0, 2*(intensity - 0.5f))
						   ,
						   (float) // GREEN
						   Math.min(1, intensity)
						   ,
						   (float) // BLUE
						   Math.min(1,
									Math.max(0, Math.max(0, 4*(intensity - 0.75f)) + 0.3*intensity)
									)
						   );
	}
	
	public static Appearance makeLandscapeAppearance() {
		return makeLandscapeAppearance(TransparencyAttributes.SCREEN_DOOR, 0.5f);
	}
	
	public static Appearance makeLandscapeAppearance(int transparencyType, float amtTransparency) {
		// amtTransparency: 0.0 = solid, 1.0 = completely transparent
		// transparencyType: TransparencyAttributes.SCREEN_DOOR, NICEST, BLENDED, NONE etc
		Appearance                     appear = new Appearance();
		PolygonAttributes          polyAttrib = new PolygonAttributes();
		TransparencyAttributes transparAttrib = new TransparencyAttributes();
		
		if (amtTransparency <= 0.0f) { // "Is the object opaque?"
			// Since the object is OPAQUE, cull back-facing polygons (default).
			// Other options are CULL_NONE and CULL_FRONT
			polyAttrib.setCullFace(PolygonAttributes.CULL_BACK);
			//polyAttrib.setCullFace(PolygonAttributes.CULL_NONE);
			transparAttrib.setTransparencyMode(TransparencyAttributes.NONE);			
		} else {
			// Don't cull any back-facing polygons if the object is
			// partially transparent
			polyAttrib.setCullFace(PolygonAttributes.CULL_NONE);
			transparAttrib.setTransparency(amtTransparency); 
			transparAttrib.setTransparencyMode(transparencyType);
		}
		appear.setPolygonAttributes(polyAttrib);
		appear.setTransparencyAttributes(transparAttrib);
		return appear;
	}
			
	/*
	void setTerrainFromGraph() {
		// Changed non-local variables: terrain, gridVerticesX, gridVerticesY
		
		CyNetwork network = Cytoscape.getCurrentNetwork();
		CyNetworkView networkView = Cytoscape.getCurrentNetworkView();
		//can't continue if either of these is null
		if (network == null || networkView == null) {
			System.err.println("Alex: ERROR: setTerrainFromGraph (LandscapeShape3D.java). Need to have a network in order to display one! Make a network and create a view, and try again.");
			randomizeTerrain();
			return;
		}
		
		float maximumX = Float.MIN_VALUE;
		float minX = Float.MAX_VALUE;
		
		float maximumY = Float.MIN_VALUE;
		float minY = Float.MAX_VALUE;
		
		int[] nodes = network.getNodeIndicesArray();
		for ( int i = 0; i < nodes.length; ++i ) {
			float x = (float) networkView.getNodeDoubleProperty( nodes[i], GraphView.NODE_X_POSITION);
			float y = (float) networkView.getNodeDoubleProperty( nodes[i], GraphView.NODE_Y_POSITION);
			
			if (x < minX) { minX = x; }
			if (x > maximumX) { maximumX = x; }
			if (y < minY) { minY = y; }
			if (y > maximumY) { maximumY = y; }
			
			// Get the minimum and maximum boundaries for the graph
		}
		
		float xSpread = maximumX - minX;
		float ySpread = maximumY - minY;
		
		float spreadAmount = 0.05f; // X amount border on every side (0.1 = 10%, 0.05 = 5%, etc)
		
		int bufferMinX = (int) Math.floor(minX - xSpread*spreadAmount);
		int bufferMinY = (int) Math.floor(minY - ySpread*spreadAmount);
		int bufferMaximumX = (int) Math.ceil(maximumX + xSpread*spreadAmount);
		int bufferMaximumY = (int) Math.ceil(maximumY + ySpread*spreadAmount);
		
		int bufferXSpread = bufferMaximumX - bufferMinX;
		int bufferYSpread = bufferMaximumY - bufferMinY;
		
		this.gridVerticesX = 250; // controls the resolution of the final image
		this.gridVerticesY = 250;
		terrain = new float[gridVerticesX][gridVerticesY];
		for (int x = 0; x < gridVerticesX; x++) {
			for (int y = 0; y < gridVerticesY; y++) {
				terrain[x][y] = 0;
			}
		}
		
		for ( int i = 0; i < nodes.length; ++i ) {
			float x = (float) networkView.getNodeDoubleProperty(nodes[i], GraphView.NODE_X_POSITION);
			float y = (float) networkView.getNodeDoubleProperty(nodes[i], GraphView.NODE_Y_POSITION);
			
			int scaledX = (int) ((x - bufferMinX) / bufferXSpread * this.gridVerticesX);
			int scaledY = (int) ((y - bufferMinY) / bufferYSpread * this.gridVerticesY);
			scaledX = (int) Math.max(Math.min((this.gridVerticesX - 1), scaledX), 0); // bound between 0 and (this.gridVerticesX-1)
			scaledY = (int) Math.max(Math.min((this.gridVerticesY - 1), scaledY), 0); // bound between 0 and (this.gridVerticesY-1)
			
			addToTerrain(scaledX, scaledY);
		}
				
	}
	
	void addToTerrain(int midX, int midY) {
		final int dist = 8; // final = "constant" in other languages. only more ominous-sounding.
		
		for (int x = midX - dist; x <= midX + dist; x++) {
			if (x < 0 || x >= this.gridVerticesX) { continue; } // skip anything out-of-bounds
			for (int y = midY - dist; y <= midY + dist; y++) {
				if (y < 0 || y >= this.gridVerticesY) { continue; } // skip anything out-of-bounds
				float distanceSq = (float) Point2D.distanceSq(x, y, midX, midY);
				terrain[x][y] += (dist * 1.5) * (dist * 1.5) - distanceSq;
				if (terrain[x][y] > maxTerrainValue) {
					maxTerrainValue = terrain[x][y]; // update the maximum terrain value if necessary
				}
			}
		}
	}
	*/
} 

