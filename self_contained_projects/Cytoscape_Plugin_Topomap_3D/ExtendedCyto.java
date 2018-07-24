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

//  ExtendedCyto.java
//  Topomap
//
//  Created by alexgraehl on 12/2/05.

// Java stuff
import javax.vecmath.*;
import java.lang.Math;
import java.util.*;
import javax.vecmath.Point3d;

// 3d things
import com.sun.j3d.utils.applet.MainFrame; 
import com.sun.j3d.utils.universe.*;
import javax.media.j3d.*;
import com.sun.j3d.utils.geometry.*; // Sphere, Cone, etc.

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

import java.awt.geom.Point2D;

public class ExtendedCyto {
	
	public static final int MIN_PT = 0;
	public static final int MAX_PT = 1;
	
	public static Point2D.Float[] calculateBounds(CyNetworkView cnv) {
		// Calculates the minimum and maximum X/Y points for the network view that is passed in
		// Returns a two-element array.
		
		float minX = Float.MAX_VALUE;
		float minY = Float.MAX_VALUE;
		
		float maximX = Float.MIN_VALUE;
		float maximY = Float.MIN_VALUE;
		
		for (Iterator it = cnv.getNodeViewsIterator(); it.hasNext(); ) {
			NodeView nView = (NodeView)it.next();
			CyNode    node = (CyNode)nView.getNode();
			
			float thisX = (float)nView.getXPosition();
			float thisY = (float)nView.getYPosition();
			
			if (thisX > maximX) { maximX = thisX; }
			if (thisY > maximY) { maximY = thisY; }
			
			if (thisX < minX) { minX = thisX; }
			if (thisY < minY) { minY = thisY; }
		}
		
		Point2D.Float[] minAndMax = new Point2D.Float[2];
		
		minAndMax[MIN_PT] = new Point2D.Float(minX, minY);
		minAndMax[MAX_PT] = new Point2D.Float(maximX, maximY);
		
		return minAndMax;
	}
	
}
