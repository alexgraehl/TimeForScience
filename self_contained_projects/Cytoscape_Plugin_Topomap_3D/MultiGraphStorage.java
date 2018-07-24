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

import java.util.Set;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;
import javax.swing.JOptionPane; // <-- for showing error/warning/info dialog boxes
import java.awt.geom.Point2D;

import giny.view.EdgeView;
import giny.view.NodeView;
import giny.view.*;

import cytoscape.Cytoscape;
import cytoscape.CyNetwork;
import cytoscape.CyNode;
import cytoscape.CyEdge;
import cytoscape.view.CyNetworkView;
import cytoscape.layout.*;
import cytoscape.*;
import cytoscape.view.CytoscapeDesktop;
import cytoscape.visual.*;
import cytoscape.visual.calculators.*;
import cytoscape.visual.mappings.*;
import cytoscape.data.Semantics;

import java.awt.Color;


public class MultiGraphStorage {
		
	private Vector cn; // CyNetwork
	private Vector cnv; // CyNetworkView
	private Vector centroids; // Point2D.Float
	private Vector nodeFillColors; // java.awt.Color: the color of nodes in the cluster at index i
	
	private int       numItems;
	private CyNetwork fullNet;
	
	public CyNetwork getOriginalNetwork() { return this.fullNet; } // the source network that started all the multi-layout-ing
	
	MultiGraphStorage(CyNetwork fullNetwork) { this(fullNetwork, 10); } // constructor, with a default value
	MultiGraphStorage(CyNetwork fullNetwork, int numGraphs) { // constructor
		// constructor
		this.numItems = 0;
		this.fullNet = fullNetwork;
		cn = new Vector(numGraphs);
		cnv = new Vector(numGraphs);
		centroids = new Vector(numGraphs);
		nodeFillColors         = new Vector(numGraphs);
	}
	
	public void setCentroid(int index, Point2D.Float cent) {
		centroids.set(index, cent);
	}
	
	public Point2D.Float getCentroid(int index) { return (Point2D.Float)centroids.get(index); }
	public CyNetwork getNetwork(int index) { return (CyNetwork) cn.get(index); }
	
	public CyNetworkView getNetworkView(int index) { return (CyNetworkView) cnv.get(index); }
	
	public java.awt.Color getNodeFillColor(int index) {	return (java.awt.Color)nodeFillColors.get(index); }
	
	public void add(Set nodeViews) {	add(nodeViews, "No name", null);	}
	public void add(Set nodeViews, String graphName) {	add(nodeViews, graphName, null);	}
	public void add(Set nodeViews, String graphName, java.awt.Color nodeColor) {
		// Adds a new graph to the graph storage, from a list of nodeViews
		Set nodeSet = new HashSet();
		Set edgeSet = new HashSet();
		
		if (nodeColor == null) {
			nodeFillColors.add(java.awt.Color.gray); // default
		} else {
			nodeFillColors.add(nodeColor);
		}
		
		// First:
		// add all the nodes into the "nodeIdxSet"
		for (Iterator it = nodeViews.iterator(); it.hasNext(); ) {
			NodeView nv = (NodeView) it.next();
			CyNode node = (CyNode) nv.getNode();
			//int nodeIndex = this.fullNet.getIndex(node);
			nodeSet.add(node); // this node needs to go in
		}
		
		// Now go through and figure out which edges we should add.
		for (Iterator it = nodeSet.iterator(); it.hasNext(); ) {
			CyNode node = (CyNode) it.next();
			java.util.List edgeList = ((giny.model.GraphPerspective)this.fullNet).getAdjacentEdgesList(node, true, true, true);
			for (Iterator edgeIt = edgeList.iterator(); edgeIt.hasNext(); ) {
				CyEdge edge = (CyEdge) edgeIt.next();
				
				CyNode sourceN = (CyNode) edge.getSource();
				CyNode targetN = (CyNode) edge.getTarget();
				
				if (nodeSet.contains(sourceN) && nodeSet.contains(targetN)) {
					edgeSet.add(edge); // this edge needs to go in ONLY if both its endpoints are in the nodeIdxSet
				}
			}
		}
		
		// Get all the IDs for the nodes
		int[] nodeIdxArr = new int[nodeSet.size()];
		int i = 0;
		for (Iterator it = nodeSet.iterator(); it.hasNext(); ) {
			CyNode n = (CyNode)it.next();
			nodeIdxArr[i++] = (int)this.fullNet.getIndex(n);
		}
		
		// Get all the IDs for the edges
		int[] edgeIdxArr = new int[edgeSet.size()];
		i = 0;
		for (Iterator it = edgeSet.iterator(); it.hasNext(); ) {
			CyEdge e = (CyEdge)it.next();
			edgeIdxArr[i++] = (int)this.fullNet.getIndex(e);
		}
		
		// createNetwork is currently BROKEN in cytoscape (I think) so that it doesn't properly handle collections of nodes.
		// That is why we have to pass in the node indices instead.
		// it is possible that createNetwork will take a Collection of objects, so this whole ID thing may be unnecessary.
		CyNetwork newNet = (CyNetwork) MultiLevelLayoutAction.createCyNetwork(nodeIdxArr, edgeIdxArr, graphName, false);
		
		int thisGroupIndex = this.size();
		for (Iterator it = newNet.nodesIterator(); it.hasNext(); ) {
			CyNode theNode = (CyNode) it.next();			
			newNet.setNodeAttributeValue(theNode, VisualStyleFactory.NODE_TYPE_ATT, VisualStyleFactory.NODE_TYPE_VALUE_REGULAR); // <-- this is a regular node
			MultiLevelLayoutAction.setGroupOfNode(thisGroupIndex, theNode); // which cluster (index) a node is in
			// Give each node an ID that tells it which SINGLE group it is a part of.
		}
		
		cn.add(newNet);
		cnv.add( (CyNetworkView)Cytoscape.createNetworkView(newNet) );
		centroids.add(null); // <-- we do this in order to keep the Vector the right size
		
		this.numItems++;
	}
	
	public int size() {
		return this.numItems;
	}
}

