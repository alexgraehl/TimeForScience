/**  Copyright (c) 2003 Institute for Systems Biology
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  The software and
**  documentation provided hereunder is on an "as is" basis, and the
**  Institute for Systems Biology has no obligations to provide maintenance, 
**  support, updates, enhancements or modifications.  In no event shall the
**  Institute for Systems Biology be liable to any party for direct, 
**  indirect, special,incidental or consequential damages, including 
**  lost profits, arising out of the use of this software and its 
**  documentation, even if the Institute for Systems Biology 
**  has been advised of the possibility of such damage. See the
**  GNU General Public License for more details.
**   
**  You should have received a copy of the GNU General Public License
**  along with this program; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**/
import cytoscape.*;
import cytoscape.view.CytoscapeDesktop;
import cytoscape.visual.*;
import cytoscape.visual.calculators.*;
import cytoscape.visual.mappings.*;
import cytoscape.data.Semantics;

import java.awt.Color;
import java.awt.Font;

public class VisualStyleFactory {
	public static final String ABSTRACT_METANODE_VS = "Meta-Node Coloring";
	
	// ============== GLOBAL constants ==============
	private static final Color   BACKGROUND_COLOR = Color.white;
	
	// ============== EDGE-related constants ==============
	public static final String EDGE_CLASS_ATT = "edgeClass";
	public static final String   EDGE_CLASS_CENTROID_TO_CENTROID = "centroidToCentroid";
	public static final String   EDGE_CLASS_WITHIN_CLUSTER = "internalCluster";
	
	public static final String   EDGE_GROUP_ATT = "edgeGroup"; // if an edge is between two nodes in the SAME group, it is part of that group, and gets assigned a group number.
	public static final String   EDGE_GROUP_CENTROID = "centroid"; // instead of a group number, centroids are in the "centroid" group
	public static final String   EDGE_GROUP_CLUSTER_TO_OUTSIDE = "clusterToOutside";

	private static final Color   DEFAULT_EDGE_COLOR = new java.awt.Color(0.7f, 0.7f, 0.7f);
	private static final Color   CENTROID_TO_CENTROID_EDGE_COLOR = Color.red;
	
	private static final LineType DEFAULT_EDGE_LINE = LineType.LINE_1;
	private static final LineType WITHIN_CLUSTER_LINE = LineType.LINE_4;
	private static final LineType CENTROID_EDGE_LINE = LineType.DASHED_5;
	
	// ====================================================
	// ============== NODE-related constants ==============
	public static final String NODE_TYPE_ATT = "nodeType";
	public static final String   NODE_TYPE_VALUE_REGULAR  = "regular";
	public static final String   NODE_TYPE_VALUE_CENTROID = "centroid";
	
	public static final String NODE_GROUP_NUM_ATT = "nodeGroup"; // which cluster (index) a node or centroid is in
	
	private static final Color DEFAULT_NODE_COLOR = Color.gray;
	
	private static final LineType DEFAULT_BORDER_LINE = LineType.LINE_1;
	private static final LineType CENTROID_BORDER_LINE = LineType.LINE_5;
	
	private static final Double DEFAULT_NODE_WIDTH    = new Double(50);
	private static final Double DEFAULT_NODE_HEIGHT   = DEFAULT_NODE_WIDTH;
	
	private static final Double CENTROID_WIDTH_VALUE  = new Double(100);
	private static final Double CENTROID_HEIGHT_VALUE = CENTROID_WIDTH_VALUE;
    private static final int NUM_NODE_COLOR_GROUPS = 10; // how many predefined color groups are there? (everything else will be the default color)
	
	private static final Integer DEFAULT_NODE_FONT_SIZE = new Integer(36);
	private static final Integer CENTROID_NODE_FONT_SIZE = new Integer(48);
	
	public static final Color DEFAULT_NODE_FONT_COLOR = Color.black;
	public static final Color CENTROID_NODE_FONT_COLOR = Color.black;
	
	// ====================================================
	
	//setNodeBorderColorCalculator(NodeColorCalculator c) 
	//void	setNodeFillColorCalculator(NodeColorCalculator c) 
	//			void	setNodeFontFaceCalculator(NodeFontFaceCalculator c) 
	//			void	setNodeFontSizeCalculator(NodeFontSizeCalculator c) 
	//			void	setNodeHeightCalculator(NodeSizeCalculator c) 
	//			void	setNodeLabelCalculator(NodeLabelCalculator c) 
	//			void	setNodeLineTypeCalculator(NodeLineTypeCalculator c) 
	//			void	setNodeShapeCalculator(NodeShapeCalculator c) 
	//			void	setNodeSizeLocked(boolean b) 
	//			void	setNodeToolTipCalculator(NodeToolTipCalculator c) 
	//			void	setNodeWidthCalculator(NodeSizeCalculator c) 
	
	public static VisualStyle createCustomVisualStyle(CyNetwork network, MultiGraphStorage mgs){
		CytoscapeDesktop cyDesktop = Cytoscape.getDesktop();
		VisualMappingManager vmManager = cyDesktop.getVizMapManager();
		NodeAppearanceCalculator nodeAppCalc = new NodeAppearanceCalculator();
		EdgeAppearanceCalculator edgeAppCalc = new EdgeAppearanceCalculator();
		CalculatorCatalog calculatorCatalog = vmManager.getCalculatorCatalog();
		GlobalAppearanceCalculator globalAppCalc = vmManager.getVisualStyle().getGlobalAppearanceCalculator();
		
		VisualStyleFactory.handleNodeCalculators(nodeAppCalc, network, mgs); // NODES
		VisualStyleFactory.handleEdgeCalculators(edgeAppCalc, network, mgs); // EDGES
		VisualStyleFactory.handleGlobalCalculators(globalAppCalc, network, mgs); // GLOBAL
		
		//------------------------- Create a visual style -------------------------------//
		VisualStyle visualStyle = new VisualStyle(ABSTRACT_METANODE_VS, nodeAppCalc, edgeAppCalc, globalAppCalc);
		//calculatorCatalog.addVisualStyle(visualStyle);
		return visualStyle;
	}
	
	// (***********************************************************************) //
	// (***********************************************************************) //
	public static void handleGlobalCalculators(GlobalAppearanceCalculator globalAppCalc, CyNetwork network, MultiGraphStorage mgs) {
		globalAppCalc.setDefaultBackgroundColor(BACKGROUND_COLOR);
	}
	
	// (***********************************************************************) //
	// (***********************************************************************) //
	public static void handleNodeCalculators(NodeAppearanceCalculator nodeAppCalc, CyNetwork network, MultiGraphStorage mgs) {
		// ------------------------------ Set the label ------------------------------//
		// Display the value for Semantics.COMMON_NAME as a label
		String cName = "Common name";
		NodeLabelCalculator nlc = (Cytoscape.getDesktop().getVizMapManager().getCalculatorCatalog()).getNodeLabelCalculator(cName);
		if (nlc == null) {
			PassThroughMapping m =
			new PassThroughMapping(new String(), Semantics.COMMON_NAME);
			nlc = new GenericNodeLabelCalculator(cName, m);
		}
		nodeAppCalc.setNodeLabelCalculator(nlc);
		
		// ------------------------------ Set node colors ---------------------------//
		DiscreteMapping colorMapping = new DiscreteMapping(DEFAULT_NODE_COLOR, ObjectMapping.NODE_MAPPING);
		colorMapping.setControllingAttributeName(NODE_GROUP_NUM_ATT, network, false);
		for (int i = 0; i < mgs.size(); i++) {
			colorMapping.putMapValue(new String("" + i), mgs.getNodeFillColor(i));
			// Fill the mappings with Color.red, Color.blue, etc.
		}
		nodeAppCalc.setNodeFillColorCalculator(new GenericNodeColorCalculator("CustomNodeColor", colorMapping));
		
		// ------------------------------ Set node shapes ---------------------------//
		Byte defaultShape = new Byte(ShapeNodeRealizer.ELLIPSE); // <-- default nodes: ellipses
		DiscreteMapping shapeMapping = new DiscreteMapping(defaultShape, ObjectMapping.NODE_MAPPING);
		shapeMapping.setControllingAttributeName(NODE_TYPE_ATT,
												 network, false);
		//shapeMapping.putMapValue("metaNode", new Byte(ShapeNodeRealizer.ELLIPSE));
		
		shapeMapping.putMapValue(NODE_TYPE_VALUE_CENTROID, new Byte(ShapeNodeRealizer.RECT)); // centroids: rectangles
		shapeMapping.putMapValue(NODE_TYPE_VALUE_REGULAR,  new Byte(ShapeNodeRealizer.ELLIPSE));
		// Possible shape values: static byte	DIAMOND ELLIPSE HEXAGON OCTAGON PARALLELOGRAM RECT RECT_3D ROUND_RECT TRAPEZOID TRAPEZOID_2 TRIANGLE 
		nodeAppCalc.setNodeShapeCalculator(new GenericNodeShapeCalculator("CustomNodeShape", shapeMapping));
		
		//---------------------- Set the thickness of the border -------------------//
		DiscreteMapping borderMapping = new DiscreteMapping(DEFAULT_BORDER_LINE, ObjectMapping.NODE_MAPPING);
		borderMapping.setControllingAttributeName(NODE_TYPE_ATT, network, false);
		borderMapping.putMapValue(NODE_TYPE_VALUE_CENTROID, CENTROID_BORDER_LINE); // centroids: thick borders
		nodeAppCalc.setNodeLineTypeCalculator(new GenericNodeLineTypeCalculator("CustomNodeBorderThickness", borderMapping));
		
		//--------------------- Set the size of the nodes --------------------------//
		DiscreteMapping wMapping = new DiscreteMapping(DEFAULT_NODE_WIDTH, ObjectMapping.NODE_MAPPING);
		wMapping.setControllingAttributeName(NODE_TYPE_ATT, network, false); // used to be "true"
		wMapping.putMapValue(NODE_TYPE_VALUE_CENTROID, CENTROID_WIDTH_VALUE);
		nodeAppCalc.setNodeWidthCalculator(new GenericNodeSizeCalculator("CustomNodeWidth", wMapping));
		
		DiscreteMapping hMapping = new DiscreteMapping(DEFAULT_NODE_HEIGHT, ObjectMapping.NODE_MAPPING);
		hMapping.setControllingAttributeName(NODE_TYPE_ATT, network, false); // used to be "true"
		hMapping.putMapValue(NODE_TYPE_VALUE_CENTROID, CENTROID_HEIGHT_VALUE);
		nodeAppCalc.setNodeHeightCalculator(new GenericNodeSizeCalculator("CustomNodeHeight", hMapping));
		
		// ------------------------------ Font sizes of labels -----------------------------//
		DiscreteMapping fontSizeMapping = new DiscreteMapping(DEFAULT_NODE_FONT_SIZE, ObjectMapping.NODE_MAPPING);
		fontSizeMapping.setControllingAttributeName(NODE_TYPE_ATT, network, false);
		fontSizeMapping.putMapValue(NODE_TYPE_VALUE_CENTROID, CENTROID_NODE_FONT_SIZE);
		nodeAppCalc.setNodeFontSizeCalculator(new GenericNodeFontSizeCalculator("CustomNodeFontSize", fontSizeMapping));
				
		// ------------------------------ Font faces of labels -----------------------------//
		DiscreteMapping fontFaceMapping = new DiscreteMapping(new Font(null, Font.PLAIN, 12), ObjectMapping.NODE_MAPPING);
		fontFaceMapping.setControllingAttributeName(NODE_TYPE_ATT, network, false);
		fontFaceMapping.putMapValue(NODE_TYPE_VALUE_CENTROID, new Font(null, Font.BOLD, 12)); // make centroid labels bold
		nodeAppCalc.setNodeFontFaceCalculator(new GenericNodeFontFaceCalculator("CustomNodeFontFace", fontFaceMapping));
	}
	
	public static void handleEdgeCalculators(EdgeAppearanceCalculator edgeAppCalc, CyNetwork network, MultiGraphStorage mgs) {
		// ------------------------------ Set edge colors ---------------------------//
		DiscreteMapping colorMapping = new DiscreteMapping(DEFAULT_EDGE_COLOR, ObjectMapping.NODE_MAPPING); // the default does NOT work
		colorMapping.setControllingAttributeName(EDGE_GROUP_ATT, network, false);
		colorMapping.putMapValue(EDGE_GROUP_CENTROID, CENTROID_TO_CENTROID_EDGE_COLOR); // centroid-centroid edges are this color
		colorMapping.putMapValue(EDGE_GROUP_CLUSTER_TO_OUTSIDE, DEFAULT_EDGE_COLOR);
		for (int i = 0; i < mgs.size(); i++) {
			colorMapping.putMapValue(new String("" + i), mgs.getNodeFillColor(i)); // within-group edges of index "i" are this color
			// Fill the mappings with Color.red, Color.blue, etc.
		}
		edgeAppCalc.setEdgeColorCalculator(new GenericEdgeColorCalculator("Set Edge Color", colorMapping));
		
		// ------------------------------ Set edge line type ---------------------------//
		DiscreteMapping lineMapping = new DiscreteMapping(DEFAULT_EDGE_LINE, ObjectMapping.NODE_MAPPING);
		lineMapping.setControllingAttributeName(EDGE_CLASS_ATT, network, false);
		lineMapping.putMapValue(EDGE_CLASS_CENTROID_TO_CENTROID, CENTROID_EDGE_LINE); // centroid-centroid edges
		lineMapping.putMapValue(EDGE_CLASS_WITHIN_CLUSTER,       WITHIN_CLUSTER_LINE); // centroid-centroid edges
		edgeAppCalc.setEdgeLineTypeCalculator(new GenericEdgeLineTypeCalculator("Set Edge Style", lineMapping));
	}
	// (***********************************************************************) //	
}