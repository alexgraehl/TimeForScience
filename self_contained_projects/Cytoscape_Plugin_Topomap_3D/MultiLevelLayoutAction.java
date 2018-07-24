import java.util.*;
import java.awt.*;
import java.awt.geom.Point2D;
import java.awt.event.ActionEvent;
import javax.swing.JOptionPane;

import cytoscape.visual.*; // <-- includes the LineType class
import cytoscape.visual.calculators.*;
import cytoscape.visual.mappings.*;
import cytoscape.visual.VisualMappingManager;
import cytoscape.visual.NodeAppearanceCalculator;

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

//import giny.model.GraphPerspective;
import giny.model.*;
//import giny.view.NodeView;
import giny.view.*;
//import giny.util.*; // the layout engines
import giny.util.Sugiyama;

import cytoscape.plugin.CytoscapePlugin;
import cytoscape.Cytoscape;
import cytoscape.CyNetwork;
import cytoscape.CyNode;
import cytoscape.CyEdge;
import cytoscape.view.CyNetworkView;
import cytoscape.data.Semantics;
//import cytoscape.data.readers.GraphReader;
import cytoscape.layout.*;
//import cytoscape.data.ExpressionData;
import cytoscape.util.CyFileFilter;
import cytoscape.util.CytoscapeAction;
import cytoscape.util.FileUtil;

import java.io.FileInputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;

import phoebe.PGraphView;

public class MultiLevelLayoutAction extends CytoscapeAction {
	
	/*public static final String nodeFillColorBypass = "node.fillColor"; 
	public static final String nodeBorderColorBypass = "node.borderColor"; 
	public static final String nodeLineTypeBypass = "node.lineType"; 
	public static final String nodeShapeBypass = "node.shape"; 
	// nodeFontBypass
	// nodeLineTypeBypass
	public static final String nodeWidthBypass = "node.width"; 
	public static final String nodeHeightBypass = "node.height"; 
	public static final String nodeLabelBypass = "node.label"; 
	public static final String nodeToolTipBypass = "node.toolTip"; 
	public static final String nodeFontBypass = "node.font"; 
	public static final String nodeLabelColorBypass = "node.labelColor";*/
		
	public static final Double BASE_NODE_SIZE = new Double(25); // normal node height / width
	public static final Double CENTROID_SIZE = new Double(BASE_NODE_SIZE.doubleValue() * 2); // height & width of centroid nodes
	public static final Float CENTROID_OPACITY = new Float(0.75f); // alpha value: 0.0 = transparent, 1.0 = opaque
	
	private static String NUM_GROUPS_KEY = new String("NUM_GROUPS_KEY");
	
	private HashMap nodeGroupMap = null;
	
	public MultiLevelLayoutAction(String menuTitle) { super(menuTitle); }
	
	//public MultiLevelLayoutAction() {super("Multi-level Layout");}
	
	// actionPerformed: this is what happens when the menu item is selected
	public void actionPerformed(ActionEvent ae) {
		System.out.println("Running a multi-level layout...");
		CyNetwork fullnet = Cytoscape.getCurrentNetwork();
		CyNetworkView fullview = Cytoscape.getCurrentNetworkView();
		if (fullnet == null || fullview == null) {
			System.out.println("MultiLevelLayout failed, due to the lack of a network and view.");
			return; // won't work if we don't have a network and view
		}
		
		String modesFilename = chooseModesFile();
		if (null == modesFilename) {
            return; // <-- if the file dialog is CANCELLED without making a selection, we get here (not an error)
		}
		
		int numGroupsFoundInModesFile = 0; // How many different MODES clusters are there?
		try	{
			nodeGroupMap = loadModesClusters(fullview, modesFilename);
			numGroupsFoundInModesFile = ((Integer)nodeGroupMap.get(MultiLevelLayoutAction.NUM_GROUPS_KEY)).intValue();
		} catch (Exception e) {
			String message = "Error in reading the MODES cluter data file \"" + modesFilename + "\".\nIt is possible that it is improperly formatted, does not have the proper permissions, or else there is some other file error.";
			JOptionPane.showMessageDialog(Cytoscape.getDesktop(), message);
			System.err.println(message);
			System.err.println(e);
			e.printStackTrace(System.err);
		}
		
		try {
			Set clusterless = new HashSet(); // <-- things without known cluster assignments
			Set[] clusters = new HashSet[numGroupsFoundInModesFile];
			for (int i = 0; i < numGroupsFoundInModesFile; i++) {
				clusters[i] = new HashSet(); // split the nodes up by group
			}
			
			for (Iterator iter = fullview.getNodeViewsIterator(); iter.hasNext(); ) {
				// Go through all the nodes...
				NodeView nView = (NodeView)iter.next();
				CyNode    node = (CyNode)nView.getNode();
				String nodeCanonicalName = (String)fullnet.getNodeAttributeValue(node, Semantics.CANONICAL_NAME);
				Integer nodeGroupID = ((Integer)nodeGroupMap.get(nodeCanonicalName)); // get the node group from the hashtable
				if (nodeGroupID == null) {
					System.out.println("Node \"" + nodeCanonicalName + "\" was assigned to \"clusterless\"...");
					clusterless.add(nView);
				} else {
					clusters[ nodeGroupID.intValue() ].add(nView);
				}
			}
			
			MultiGraphStorage mg = new MultiGraphStorage(fullnet);
			for (int i = 0; i < clusters.length; i++) {
				mg.add(clusters[i], "G" + i, getColorForGraph(i));
			}
			mg.add(clusterless, "Clusterless", Color.gray); // clusterless nodes are gray
			
			CyNetworkView newNetView = MultiLevelLayoutAction.makeClusterNetworkView(mg);
			CyNetwork     newNet = newNetView.getNetwork();
			//Cytoscape.setCurrentNetwork(newNet.getIdentifier()); // <-- deprecated
			
			VisualMappingManager vMapManager = Cytoscape.getDesktop().getVizMapManager(); 
			vMapManager.setVisualStyle( (cytoscape.visual.VisualStyle)VisualStyleFactory.createCustomVisualStyle(newNet, mg) );
			newNetView.redrawGraph(true, true);
			newNetView.fitContent();
		} catch (Exception e) {
			String message = "An error occurred while running the layout.";
			JOptionPane.showMessageDialog(Cytoscape.getDesktop(), message);
			System.err.println(e);
			e.printStackTrace(System.err);
		}
	} // void actionPerformed()
	
	public static void setGroupOfEdge(int groupIndex, CyEdge whichEdge) {
		Cytoscape.setEdgeAttributeValue(whichEdge, VisualStyleFactory.EDGE_GROUP_ATT, new String("" + groupIndex));
	}
	public static void setGroupOfEdge(int groupIndex, int whichEdgeIndex) {
		Cytoscape.setEdgeAttributeValue(Cytoscape.getRootGraph().getEdge(whichEdgeIndex), VisualStyleFactory.EDGE_GROUP_ATT, new String("" + groupIndex));
	}
	public static int getGroupOfEdge(CyEdge whichEdge) {
		return Integer.parseInt((String)Cytoscape.getEdgeAttributeValue(whichEdge, VisualStyleFactory.EDGE_GROUP_ATT));
	}
	public static int getGroupOfEdge(int edgeIndex) {
		return Integer.parseInt((String)Cytoscape.getEdgeAttributeValue(Cytoscape.getRootGraph().getEdge(edgeIndex), VisualStyleFactory.EDGE_GROUP_ATT));
	}
	
	public static void setGroupOfNode(int groupIndex, Node whichNode) {
		Cytoscape.setNodeAttributeValue(whichNode, VisualStyleFactory.NODE_GROUP_NUM_ATT, new String("" + groupIndex));
	}
	public static void setGroupOfNode(int groupIndex, int whichNodeIndex) {
		Cytoscape.setNodeAttributeValue(Cytoscape.getRootGraph().getNode(whichNodeIndex), VisualStyleFactory.NODE_GROUP_NUM_ATT, new String("" + groupIndex));
	}
	public static int getGroupOfNode(Node whichNode) {
		return Integer.parseInt((String)Cytoscape.getNodeAttributeValue(whichNode, VisualStyleFactory.NODE_GROUP_NUM_ATT));
	}
	public static int getGroupOfNode(int nodeIndex) {
		return Integer.parseInt((String)Cytoscape.getNodeAttributeValue(Cytoscape.getRootGraph().getNode(nodeIndex), VisualStyleFactory.NODE_GROUP_NUM_ATT));
	}
	
	private static String chooseModesFile() {
		// Chooses a file with Modes (cluster) data and returns the filename
		String returnedFilename;
		try {
			CyFileFilter filter = new CyFileFilter(); // <-- default filter accepts ALL file types
			filter.setDescription("MODES cluster output files");
            returnedFilename = FileUtil.getFile("Load MODES cluster groupings file",
												FileUtil.LOAD,
												new CyFileFilter[]{filter}).toString();
        } catch (Exception e) {
			System.out.println("Multi-level layout was cancelled: no MODES file was selected...");
			returnedFilename = null;
        }
		return returnedFilename;
	}
	
	public static float fractPart(float x) { return ((float)x - (float)Math.floor(x)); } // returns the non-integer part of a float
	public static float fractPart(double x) { return ((float)x - (float)Math.floor(x)); } // returns the non-integer part of a double
	
	public static Color getColorForGraph(int index) {
		if (index < 10) {
			float hue        = fractPart(index*0.3);
			float saturation = 0.9f;
			float bright     = 0.9f;
			return Color.getHSBColor(hue, saturation, bright);
		} else if (index < 20) {
			float hue        = fractPart(index*0.3 + 0.15);
			float saturation = 0.7f;
			float bright     = 0.6f;
			return Color.getHSBColor(hue, saturation, bright);
		} else {
			return Color.gray;
		}
	}
	
	public static HashMap loadModesClusters(CyNetworkView netView, String modesFilename) throws Exception {
		// 1. Format of a MODES cluster file:
		//    GroupID  NumberGenes  NumEdges  List of genes
		//    0  <tab> 4    <tab>   5  <tab>  alpha-1 <tab> alpha-2 <tab> alpha-3 <tab> alpha-4   <newline>
		//    1        2            1         beta-1 <tab> beta-2   <newline>
		//    2        1            0         gamma-1        <newline>
		
		// 2. We can ALSO accept a file of the following format:
		//    GroupID  List of genes
		//    a  <tab> alpha-1   <tab>   alpha-2   <newline>
		//    b        beta-1    <newline>
		//    c        gamma-1   <newline>
		
		final int NUM_SPECIAL_FIELDS        = 3; // all the fields except for the gene name
		final int GROUP_NAME_INDEX          = 0;
		final int NUM_GENES_IN_GROUP_INDEX  = 1;
		final int NUM_EDGES_IN_GROUP_INDEX  = 2;
		
		int numberOfGroups = 0;
		
		CyNetwork theNet = netView.getNetwork();
		
		int expectedNumberOfNodes = netView.getNetwork().getNodeCount();
		HashMap nodeGroupAssignments = new HashMap(expectedNumberOfNodes + 1);
		
		FileInputStream fstream = new FileInputStream(modesFilename);
		
		//assert (h != null) : "h was null in setTerrainFromDensity";
		BufferedReader data = new BufferedReader(new InputStreamReader(fstream));
		for (String theLine = data.readLine(); theLine != null; theLine = data.readLine()) { // read each line...
			String[] modesTokens = theLine.split("\t");
			
			if (modesTokens.length < 1) {
				continue; // probably it's a blank line
			}
			
			String groupName = modesTokens[GROUP_NAME_INDEX];
			int numGenes, numEdges;
			
			// There are two options now: the line might contain three fields at the beginning (groupName, numGenes, numEdges),
			// or it might just contain (groupName)
			// In either case, the items afterward will all be gene names.
			
			//System.out.println(theLine);
			
			int startReadingGenesAtIndex = 3; // <-- default: the 4th item (index 3) is where the gene names start, in a normal MODES file. However, sometimes the file just consists of a groupName and the gene list, in which case the start index will be 1. That's what we check below...
			try {
				numGenes = Integer.parseInt(modesTokens[NUM_GENES_IN_GROUP_INDEX]);
				numEdges = Integer.parseInt(modesTokens[NUM_EDGES_IN_GROUP_INDEX]);
			} catch (Exception e) {
				// If this FAILS, then it means we didn't have two integers at indices 1 and 2 at the beginning,
				// so we just have (groupName) (gene1) (gene2) (gene3).
				startReadingGenesAtIndex = 1;
			}
			
			for (int i = startReadingGenesAtIndex; i < modesTokens.length; i++) {
				String nameOfGene = modesTokens[i];
				if (nodeGroupAssignments.get(nameOfGene) != null) {
					System.out.println("Multiple cluster assignments for: " + nameOfGene);
				}
				
				nodeGroupAssignments.put(nameOfGene, new Integer(numberOfGroups)); // pairing is (geneName / groupIndex)
				System.out.println("Adding " + nameOfGene + " to group \"" + numberOfGroups + "\"" + " (group name in file was \"" + groupName + "\")");
			}
			numberOfGroups++;
		}
		data.close();
		
		if (numberOfGroups < 1) {
			throw new Exception("We didn't end up reading any groups from the MODES file. Something went wrong.");
		}
		
		nodeGroupAssignments.put(MultiLevelLayoutAction.NUM_GROUPS_KEY, new Integer(numberOfGroups));
		
		return nodeGroupAssignments;
	}
		
	private static CyNode createNewCentroidNode(int whichGroupIndex) {
		return Cytoscape.getCyNode("Centroid " + whichGroupIndex, true);
	}
	
	private static CyNode getCentroid(int whichGroupIndex) {
		return Cytoscape.getCyNode("Centroid " + whichGroupIndex);
	}
	
	private static float calculateMaxRadius(Point2D.Float centerPt, Iterator nvIter) {
		// Figures out the farthest point from the centerPt, and returns that distance.
		double maxRadiusSquared = 0;
		while (nvIter.hasNext()) {
			// Go through all the nodes...
			NodeView nView = (NodeView)nvIter.next();
			//System.out.println("Position of original node: " + (float)nView.getYPosition() + ", " + (float)nView.getXPosition() + ".");
			double distanceSqToThisPoint = centerPt.distanceSq(new Point2D.Float((float)nView.getXPosition(), (float)nView.getYPosition()));
			if (distanceSqToThisPoint > maxRadiusSquared) { maxRadiusSquared = distanceSqToThisPoint; }
		}
		return ( (float) Math.sqrt(maxRadiusSquared) );
	}
	
	public static void scaleNodesToRadius(float radius, Point2D.Float centerPt, CyNetworkView netView) {
		// Scales the nodes so that they fit exactly within "desiredRadius" of the origin.
		Point2D.Float thePt = new Point2D.Float(0.0f, 0.0f);
		float scaleFactor = scaleFactorForDesiredRadius(radius, centerPt, netView.getNodeViewsIterator());
		clusterAroundPoint(thePt, thePt, scaleFactor, netView.getNodeViewsIterator(), netView);
	}
	
	public static void scaleNodesByFactor(float scaleFactor, CyNetworkView netView) {
		// Scales the nodes so that they fit exactly within "desiredRadius" of the origin.
		Point2D.Float thePt = new Point2D.Float(0.0f, 0.0f);
		clusterAroundPoint(thePt, thePt, scaleFactor, netView.getNodeViewsIterator(), netView);
	}
	
	private static void setNodeProperties(NodeView nv, Integer shape, Float opacity, Double width, Double height,
										  Color unselectedFillColor, Color selectedFillColor, LineType borderLineType, Color borderColor, String toolTip) {
		// Sets node visual properties.
		// Note: borderLineType can be DASHED_1 ... _5, DOTTED_1 .. _5, or LINE_1 ... _7.
		
		CyNode node = (CyNode)nv.getNode();
		if (opacity != null) nv.setTransparency(opacity.floatValue()); // can be set directly
		if (width != null)   Cytoscape.setNodeAttributeValue(node, NodeAppearanceCalculator.nodeWidthBypass, width); // setWidth
		if (height != null)  Cytoscape.setNodeAttributeValue(node, NodeAppearanceCalculator.nodeHeightBypass, height); // setHeight
		if (unselectedFillColor != null) Cytoscape.setNodeAttributeValue(node, NodeAppearanceCalculator.nodeFillColorBypass, unselectedFillColor); // setUnselected
		if (selectedFillColor != null)   nv.setSelectedPaint(selectedFillColor); // can be set directly
		if (borderColor != null)         Cytoscape.setNodeAttributeValue(node, NodeAppearanceCalculator.nodeBorderColorBypass, borderColor);
		if (borderLineType != null)      Cytoscape.setNodeAttributeValue(node, NodeAppearanceCalculator.nodeLineTypeBypass, borderLineType);
		
		//Cytoscape.setNodeAttributeValue(node, NodeAppearanceCalculator.nodeLabelColorBypass, Color.black);
		
		if (toolTip != null)         nv.setToolTip(toolTip);
	}
	
	private static void runSpringLayout(CyNetworkView theNetView) {
		//AttributeLayout alay = new AttributeLayout(theNetView);
		//alay.doLayout("");
		CustomSpringLayouter theLayout = new CustomSpringLayouter(theNetView);
		theLayout.doLayout();
	}
	
	public static void recenterAroundCurrentCentroids(CyNetworkView theNetView) {
		// Recenters the points around their corresponding "centroid" after the centroids are moved (the centroid may not actually be a centroid anymore, though, but the corresponding points will still move to cluster around it).
		
		CyNetwork  theNet;
		int centroidCount = 0;
		try {
			theNet = theNetView.getNetwork();
			while (MultiLevelLayoutAction.getCentroid(centroidCount) != null) {
				centroidCount++;
			}
			if (centroidCount <= 0) { throw new Exception("There were no centroids on the graph, so we can't do any recentering."); }
		} catch (Exception e) {
			System.err.println(e);
			JOptionPane.showMessageDialog(Cytoscape.getDesktop(), e.getMessage());
			return;
		}
		
		Vector regularNodeGroups = new Vector(centroidCount);
		for (int i = 0; i < centroidCount; i++) {
			regularNodeGroups.add(new Vector());
		}
		
		for (Iterator it = theNetView.getNodeViewsIterator(); it.hasNext(); ) {
			NodeView nView = (NodeView)it.next();
			CyNode node    = (CyNode)nView.getNode();
			if (isRegularNode(node)) { // is this a REGULAR node (i.e., not a centroid)?
				int nodeGroupIdx = getGroupOfNode(node);
				((Vector)regularNodeGroups.get(nodeGroupIdx)).add(nView);
			}
		}
		
		for (int i = 0; i < centroidCount; i++) {
			NodeView theCentroidView = theNetView.getNodeView( getCentroid(i) );
			Point2D.Float sourceCenter = calculateCentroidPt(((Vector)regularNodeGroups.get(i)).iterator());
			Point2D.Float newCenter = new Point2D.Float((float)theCentroidView.getXPosition(), (float)theCentroidView.getYPosition());
			clusterAroundPoint(sourceCenter,
							   newCenter,
							   1.0f,
							   ((Vector)regularNodeGroups.get(i)).iterator(),
							   theNetView
							   );
		}
		
		selectAnyGroupsWithSelectedCentroids(theNetView); // select the nodes so the user can move them if desired
		
	}
	
	public static boolean isRegularNode(CyNode node) {
		String nodeType  = (String)Cytoscape.getNodeAttributeValue(node, VisualStyleFactory.NODE_TYPE_ATT);
		return (VisualStyleFactory.NODE_TYPE_VALUE_REGULAR.equals(nodeType));
	}
	
	public static boolean isCentroidNode(CyNode node) {
		String nodeType  = (String)Cytoscape.getNodeAttributeValue(node, VisualStyleFactory.NODE_TYPE_ATT);
		return (VisualStyleFactory.NODE_TYPE_VALUE_CENTROID.equals(nodeType));
	}
	
	public static void selectAnyGroupsWithSelectedCentroids(CyNetworkView theNetView) {
		// Checks through all the selected nodes. For each selected centroid, we ALSO select the nodes in its same cluster.
		// Does not UNselect any nodes.
		CyNetwork theNet = theNetView.getNetwork();
		for (Iterator it = theNetView.getSelectedNodes().iterator(); it.hasNext(); ) {
			NodeView nv = (NodeView)it.next();
			CyNode node = (CyNode)nv.getNode();
			if (isCentroidNode(node)) {
				int centroidGroupIndex = getGroupOfNode(node);
				System.out.println("Selecting groups of index " + centroidGroupIndex);
				selectNodesOfGroup(centroidGroupIndex, theNetView);
			}
		}
	}
	
	public static void selectNodesOfGroup(int groupToSelectIndex, CyNetworkView theNetView) {
		// Selects all the nodes (including centroids) of a given group.
		if (MultiLevelLayoutAction.getCentroid(groupToSelectIndex) == null) {
			String warningMessage = "We tried to select a group of index " + groupToSelectIndex + ", but that index was out of bounds.";
			System.err.println(warningMessage);
			JOptionPane.showMessageDialog(Cytoscape.getDesktop(), warningMessage);
			return;
		}
		CyNetwork theNet = theNetView.getNetwork();
		for (Iterator it = theNetView.getNodeViewsIterator(); it.hasNext(); ) {
			NodeView nv = (NodeView)it.next();
			CyNode node = (CyNode) nv.getNode();
			if (getGroupOfNode(node) == groupToSelectIndex) {
				//nv.setSelected(true); // NodeViews must be individually selected based on their group.
				theNet.setFlagged(node, true);
			}
		}
		//theNetView.setSelected(selectUs); <-- does not work as expected: nodes APPEAR selected, but do not move in a group
		//Cytoscape.setFlaggedNodes(nodesToSelect, true);
		theNetView.redrawGraph(false, true);
	}
	
	private static float scaleFactorForDesiredRadius(float desiredRadius, Point2D.Float originalCenter, Iterator nodeViewsIterator) {
		float maxRadius = calculateMaxRadius(originalCenter, nodeViewsIterator);
		float positionScaleFactor = (float)desiredRadius / (float)maxRadius; // set the scale factor so that the farthest-out data point will now be "desiredRadius" distance units from the destCenter.
		return positionScaleFactor;
	}
	
	public static void clusterAroundPoint(Point2D.Float sourceCenter, Point2D.Float destCenter,
										  float positionScaleFactor, Iterator sourceNVIterator, CyNetworkView destView) {		
		// Takes all the nodes in "sourceView" and clusters them around "destCenter" based on their relative positions around "sourceCenter." Doesn't add any nodes, just moves them around.
		// The maximum radius of the farthest-away-from-the-center point is specified by "destCenter."
		
		// radius is the radius from the center point to the farthest-out point
		// sourceCenter is the old center of the graph
		// destCenter is the new center where we will move all the points to
		// sourceCenter and destCenter can be the same point, but in this case the only thing that will have any effect is the scaling factor from the "radius" term
		
		while (sourceNVIterator.hasNext()) {
			// Go through all the nodes...
			NodeView sourceNodeView = (NodeView)sourceNVIterator.next();
			CyNode   sourceNode     = (CyNode) sourceNodeView.getNode();
			float relX = (float)sourceNodeView.getXPosition() - sourceCenter.x; // relative X position, compared to the center
			float relY = (float)sourceNodeView.getYPosition() - sourceCenter.y; // relative Y position, compared to the center			
			
			NodeView destNodeView = (NodeView)destView.getNodeView(sourceNode); // Get the NodeView associated with the DESTINATION view, for the node from the SOURCE
			
			destNodeView.setXPosition(destCenter.x + (relX * positionScaleFactor)); // position CAN be set directly
			destNodeView.setYPosition(destCenter.y + (relY * positionScaleFactor));
			
			//System.out.println("sourceCenter: " + sourceCenter.x + ", " + sourceCenter.y + " and destCenter: " + destCenter.x + ", " + destCenter.y + ". PositionScaleFactor was " + positionScaleFactor + " and maxRadius was " + maxRadius + ".");
			//System.out.println("Setting node position to: (" + (destCenter.x + (relX * positionScaleFactor)) + ", " + (destCenter.y + (relY * positionScaleFactor)) + ").");
		}
	}
	
	private static java.awt.geom.Point2D.Float calculateCentroidPt(Iterator theNodeViewsIterator) {
		// Calculates the centroid of all the nodes in this network view.
		// It really only makes sense to call this AFTER the network has been laid out.
		float   xSum = 0;
		float   ySum = 0;
		int numNodes = 0;
		while (theNodeViewsIterator.hasNext()) {
			NodeView nView = (NodeView)theNodeViewsIterator.next();
			xSum += (float) nView.getXPosition();
			ySum += (float) nView.getYPosition();
			numNodes++;
		}
		return new java.awt.geom.Point2D.Float((float)xSum/numNodes, (float) ySum/numNodes);
	}
	
	private static void addEdgesBetweenCentroids(CyNetwork sourceNet, CyNetworkView centNetView) {
		// nodeIndices: the indices for the centroid nodes that are to be joined
		// sourceNet: the source network where we get the weight data
		// centNetView: the network & view that we're going to add the edges to
		int centroidCount = 0;
		while (MultiLevelLayoutAction.getCentroid(centroidCount) != null) {
			centroidCount++;
		}
		
		float[][] edgeWeights = new float[centroidCount][centroidCount];
		
		CyNetwork centNetwork = centNetView.getNetwork();
		
		// Tally ALL the weights. Yes, it counts all edges twice.
		for (Iterator it = sourceNet.nodesIterator(); it.hasNext(); ) {
			CyNode node = (CyNode) it.next();
			java.util.List edgeList = ((giny.model.GraphPerspective)sourceNet).getAdjacentEdgesList(node, true, true, true);
			for (Iterator edgeIt = edgeList.iterator(); edgeIt.hasNext(); ) {
				CyEdge edge = (CyEdge) edgeIt.next();
				CyNode sourceN = (CyNode) edge.getSource();
				CyNode targetN = (CyNode) edge.getTarget();
				
				try {
					int sourceGroupIdx, destGroupIdx;
					try {
						sourceGroupIdx = getGroupOfNode(sourceN);
						destGroupIdx   = getGroupOfNode(targetN);
					} catch (Exception e) {
						System.out.println("Error in attempting to get a node attribute value!");
						continue;
					}
					float thisEdgeWeight = 1.0f; // set the edge weights between centroids
					edgeWeights[sourceGroupIdx][destGroupIdx] += (float)thisEdgeWeight;
					//edgeWeights[destGroupIdx][sourceGroupIdx] = edgeWeights[sourceGroupIdx][destGroupIdx]; // <-- we could probably take this out and verify that the two indices [i][j] and [j][i] ended up the same (they should if we go through each node!), but right now we do not trust that.
				} catch (Exception e) {
					String warningMessage = new String("We encountered an error in the procedure for adding edges between centroids.\nOne of the indices that was checked for was not valid.\n(Either the sourceGroupIdx or destGroupIdx was invalid.)");
					System.out.println(warningMessage + "\n" + e.getMessage());
					JOptionPane.showMessageDialog(Cytoscape.getDesktop(), warningMessage);
				}
			}
		}
		
		for (int i = 0; i < centroidCount; i++) {
			for (int j = 0; j < i; j++) {
				// Add an undirected edge between each pair of nodes
				CyEdge centEdge = Cytoscape.getCyEdge(getCentroid(i), getCentroid(j), Semantics.INTERACTION, new String("pp"), true);
				centNetwork.addEdge(centEdge);
								
				EdgeView theEdgeView = centNetView.getEdgeView(centEdge);
				
				// SET VISUAL PROPERTIES FOR CENTROID-to-CENTROID EDGES
				theEdgeView.setToolTip("Centroid " + i + " to " + j + " weights: outgoing: " + edgeWeights[i][j] + ", incoming: " + edgeWeights[j][i]);
				centNetwork.setEdgeAttributeValue(centEdge, VisualStyleFactory.EDGE_GROUP_ATT, VisualStyleFactory.EDGE_GROUP_CENTROID); // <-- this is a CENTROID-to-CENTROID edge
				centNetwork.setEdgeAttributeValue(centEdge, VisualStyleFactory.EDGE_CLASS_ATT, VisualStyleFactory.EDGE_CLASS_CENTROID_TO_CENTROID);
				//theEdgeView.setSelectedPaint(Color.red);
				//theEdgeView.setUnselectedPaint(Color.orange);
				// DONE SETTING PROPERTIES FOR CENTROID-to-CENTROID EDGES
			}
		}
	}
	
	public static GraphPerspective createCyNetwork(String netName, boolean displayInGUI) {
		return createCyNetwork(null, null, netName, displayInGUI);
	}
	public static GraphPerspective createCyNetwork(int[] nodeIdxArr, int[] edgeIdxArr, String netName, boolean displayInGUI) {
		// Rowan Christmas: a "Network" is a "giny.model.GraphPerspective" with some additional methods that tie it to the DataModel. 
		// if you want to create an intermediate, use Cytoscape.getRootGraph().createGraphPerspective() (and related methods). 
		// This would be a way to make a network WITHOUT having it show up as a selectable item.
		if (displayInGUI) {
			return Cytoscape.createNetwork(nodeIdxArr, edgeIdxArr, netName);
		} else {
			giny.model.Node[] nullNode = null;	
			// Rowan Christmas: a "Network" is a "giny.model.GraphPerspective" with some additional methods that tie it to the DataModel.
			// if you want to create an intermediate, use Cytoscape.getRootGraph().createGraphPerspective() (and related methods). 
			// This would be a way to make a network WITHOUT having it show up as a selectable item.
			if (displayInGUI) {
				return Cytoscape.createNetwork(netName);
			} else {
				CyNetwork newNet = ((Cytoscape.getRootGraph()).createNetwork(nodeIdxArr, edgeIdxArr));
				newNet.setTitle(netName);
				return newNet;
			}
		}
	}
	
	// ==================================================================
	public static CyNetworkView makeClusterNetworkView(MultiGraphStorage mgs) {
		CyNetwork clusterNet = (CyNetwork) createCyNetwork("MetaNet", false);
		
		float overallNetworkRadius = 600.0f;
		//CyNode centroids[] = new CyNode[mgs.size()];
		
		for (int groupNum = 0; groupNum < mgs.size(); groupNum++) {
			// FIRST: LAYOUT EACH CLUSTER (so we can figure out where the centroid will be)
			// Spring-layout the sub-group (cluster) index index.		
			MultiLevelLayoutAction.runSpringLayout(mgs.getNetworkView(groupNum));
			
			// NOW: CALCULATE & ADD THE CENTROID of the cluster
			mgs.setCentroid(groupNum, calculateCentroidPt(mgs.getNetworkView(groupNum).getNodeViewsIterator())); // calculate and set the centroid of the graph
			
			CyNode newCentroid = createNewCentroidNode(groupNum); // true: "create a new node"
			//centroids[groupNum] = newCentroid; //Cytoscape.getRootGraph().createNode(); // Make the centroid node
			clusterNet.addNode(newCentroid); // Add the centroid to the metaNet
			
			Cytoscape.setNodeAttributeValue(newCentroid, VisualStyleFactory.NODE_TYPE_ATT, VisualStyleFactory.NODE_TYPE_VALUE_CENTROID); // inform the centroid that it is in fact a centroid
			MultiLevelLayoutAction.setGroupOfNode(groupNum, newCentroid); // set the centroid's group
			Cytoscape.setNodeAttributeValue(newCentroid, NodeAppearanceCalculator.nodeLabelColorBypass, Color.red);
		} // Now the centroids exist, but not yet laid out with respect to one another.
		
		CyNetworkView clusterNetworkView = Cytoscape.createNetworkView(clusterNet); // new network view...
		
		addEdgesBetweenCentroids(mgs.getOriginalNetwork(), clusterNetworkView); // ADD EDGES BETWEEN ALL CENTROIDS
		
		MultiLevelLayoutAction.runSpringLayout(clusterNetworkView); // <-- lay out all the centroids
		
		MultiLevelLayoutAction.scaleNodesToRadius(overallNetworkRadius, calculateCentroidPt(clusterNetworkView.getNodeViewsIterator()), clusterNetworkView); // just scales the CENTROIDS, since that's all that's on the network right now.
		
		// Now the centroids and their edges are on the graph, and laid out, but no other points are there.
		
		// Now we need to add all the other data points that are in each cluster. 
		for (int i = 0; i < mgs.size(); i++) {
			System.out.println("Running clustering loop with i = " + i);
			for (Iterator iter = mgs.getNetwork(i).nodesIterator(); iter.hasNext(); ) {
				CyNode node = (CyNode) iter.next();
				clusterNet.addNode(node); // add all the nodes from the other networks to this network
				NodeView clusterNV = (NodeView)clusterNetworkView.getNodeView(node);
				Cytoscape.setNodeAttributeValue(node, NodeAppearanceCalculator.nodeLabelColorBypass, VisualStyleFactory.DEFAULT_NODE_FONT_COLOR);
				//MultiLevelLayoutAction.setNodeProperties(clusterNV, null, null, null, null, null, null, null, null, new String("Group " + i));
			}
			NodeView thisCentroidView = clusterNetworkView.getNodeView( getCentroid(i) );
			
			Point2D.Float newCenterPt = new Point2D.Float((float)thisCentroidView.getXPosition(), (float)thisCentroidView.getYPosition());
			float radius = 1.00f * (float) Math.sqrt(mgs.getNetwork(i).getNodeCount()) * (overallNetworkRadius / mgs.size()); // controls how far apart the cluster's member nodes will appear on the graph (bigger value for radius -> nodes are more spread out)
			// Now that we added the cluster at index i's nodes, we should lay them out nicely around their corresponding centroid
			float scaleFactor = scaleFactorForDesiredRadius(radius, mgs.getCentroid(i), mgs.getNetworkView(i).getNodeViewsIterator());
			MultiLevelLayoutAction.clusterAroundPoint(mgs.getCentroid(i), // old centroid
													  newCenterPt, // new point should be the centroid location in the new graph
													  scaleFactor,
													  mgs.getNetworkView(i).getNodeViewsIterator(), // Source view
													  clusterNetworkView);    // Destination view
		}
		
		for (Iterator iter = mgs.getOriginalNetwork().edgesIterator(); iter.hasNext(); ) {
			CyEdge edge = (CyEdge)iter.next(); // Include all the edges that were present in the original graph
			clusterNet.addEdge(edge);
			try {
				CyNode sourceN = (CyNode)edge.getSource();
				CyNode targetN = (CyNode)edge.getTarget();
				int sourceGroupIdx = getGroupOfNode(sourceN);
				int destGroupIdx   = getGroupOfNode(targetN);
				if (sourceGroupIdx == destGroupIdx) {
					// source and destination clusters are the same for this edge, so let's assign it to this group
					setGroupOfEdge(sourceGroupIdx, edge);
					clusterNet.setEdgeAttributeValue(edge, VisualStyleFactory.EDGE_CLASS_ATT, VisualStyleFactory.EDGE_CLASS_WITHIN_CLUSTER);
				} else {
					// source and destination clusters are the DIFFERENT for this edge, so don't assign it to any particular group
					clusterNet.setEdgeAttributeValue(edge, VisualStyleFactory.EDGE_GROUP_ATT, VisualStyleFactory.EDGE_GROUP_CLUSTER_TO_OUTSIDE); // assign it to the "cluster to outside" group
				}
			} catch (Exception e) {
				System.out.println("Error in attempting to get an edge attribute value!");
				continue;
			}
		}
		
		return (CyNetworkView)clusterNetworkView;
	}
	// ==================================================================
}

// -----=====-----====-----=-=-=-=-=

//giny.util.Sugiyama su = new Sugiyama(view);

//su.layout();
//view = (CyNetworkView) springLayout.getGraphView();

// = new giny.util.JUNGSpringLayout(view);
//view.applyLayout((LayoutAlgorithm)new cytoscape.layout.SpringEmbeddedLayouter(view));

//view.redrawGraph(true, false);
//view.redrawGraph(false, false);


//giny.util.SpringEmbeddedLayouter springLayout = new SpringEmbeddedLayouter(view);
//cytoscape.layout.SpringEmbeddedLayouter layout = new cytoscape.layout.SpringEmbeddedLayouter(view);
//CustomSpringLayouter layout = new CustomSpringLayouter(view);
//layout.doLayout();


// -=-=-=-=-=-=-=-=-=-=-=-=-
// from google groups on cytoscape:
//Yes, you would need two node types. It could look like this: 
//DiscreteMapping disMapping = new DiscreteMapping(new Byte( 
//														   ShapeNodeRealizer.RECT),ObjectMapping.NODE_MAPPING); 
//disMapping.setControllingAttributeName("nodeType", 
//									   network, 
//									   false); 
//disMapping.putMapValue("pathway", new Byte(ShapeNodeRealizer.ELLIPSE)); 
//disMapping.putMapValue("reaction", new Byte(ShapeNodeRealizer.RECTANGLE)); 
//GenericNodeShapeCalculator shapeCalculator = 
//new GenericNodeShapeCalculator("Pathway/Reaction Shape Calculator", 
//							   disMapping); 
//nodeAppCalc.setNodeShapeCalculator(shapeCalculator); 
//
//Because in your class it seems as though you are defining the 
//
//> properties for the meta-nodes only. 
//> Secondly, is the idea I'm having now correct? 
//> My idea is, that I should get the VisualMappingManager for the network. 
//> And then create a NodeAppearanceCalculator and EdgeAppearanceCalculator 
//> with all the values needed/wanted. I cannot just change what I wanted 
//> to change and let the rest be default can I? 
//
//
//Yes, you can modify the NodeAppearanceCalculator/EdgeAppearanceCalculator 
//that is set in the current visual style: 
//VisualMappingManager manager = Cytoscape.getDesktop().getVizMapManager(); 
//NodeAppearanceCalculator nac = manager.getVisualStyle 
//().getNodeAppearanceCalculator(); 
//


/*
 //now go through our container and select each view
 for (Iterator i = nvSet.iterator(); i.hasNext(); ) {
	 NodeView nv = (NodeView)i.next();
	 //nv.setSelected(true);
	 //nv.setXPosition(55);
	 //nv.setYPosition(25);
	 //nv.setTransparency(0.5f);
	 nv.setHeight(55);
	 nv.setWidth(25);
	 nv.setSelectedPaint(Color.black);
	 nv.setUnselectedPaint(Color.yellow);
	 nv.setBorderPaint(Color.blue);
	 nv.setBorderWidth(5.0f);
 }
 */