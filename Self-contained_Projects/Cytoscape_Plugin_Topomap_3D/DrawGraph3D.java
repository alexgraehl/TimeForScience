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


// windowing stuff
import java.awt.event.*;
import java.awt.GraphicsConfiguration;
import javax.swing.JOptionPane;

// Java stuff
import javax.vecmath.*;
import java.lang.Math;
import java.util.*;
import javax.vecmath.Point3d;

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

import java.awt.geom.Point2D;

import java.util.HashMap;

public class DrawGraph3D {
	
	private static final Color3f nodeDefaultColor = new Color3f(1.0f, 1.0f, 1.0f);
	
	public static ReattachableBranchGroup newGraph(CyNetworkView theNetView) {
		return newGraph(theNetView, 1.0);
	}
	
	public static ReattachableBranchGroup newGraph(CyNetworkView theNetView, double sizeOnScreen) {
		return newGraph(theNetView, sizeOnScreen, null);
	}
	
	public static ReattachableBranchGroup newGraph(CyNetworkView theNetView, double sizeOnScreen, HashMap rememberBranches) {
		// Makes a new graph object for drawing the graph (nodes & edges, not the topographical 3d view) on the screen
		
		// sizeOnScreen: tells us how large the object should be drawn on the screen.
		// It would be 1.0 if we want it to be the exact same size as the topomap,
		// but the topomap also includes some BORDER squares which do not use
		// any graph data. If there are any border squares in the topomap, then we
		// will want to draw the graph in a slightly smaller area in order to properly
		// overlap with the topomap.
		
		//System.err.println("Size on screen for new graph was " + sizeOnScreen);
		
		ReattachableBranchGroup graphGroup = new ReattachableBranchGroup();
		// Heirarchy:
		// graphGroup_
		//   |        \_
		// nodeGroup    edgeGroup
		//
		
		if (theNetView == null) {
			System.err.println("Warning: DrawGraph3D.java: we tried to call a newGraph with a null argument.\n");
			return graphGroup;
		}
		
		Point2D.Float[] sourceGraphBounds = ExtendedCyto.calculateBounds(theNetView);
		
		ReattachableBranchGroup nodeRBGroup 
			= createNodeGroup(theNetView, sourceGraphBounds, sizeOnScreen, rememberBranches);
		
		ReattachableBranchGroup edgeRBGroup 
			= createEdgeGroup(theNetView, sourceGraphBounds, sizeOnScreen);
		if (rememberBranches != null) {
			rememberBranches.put(Topomap3D.kGRAPH_EDGES, edgeRBGroup);
		}
		
		addLightsToBranchGroup(graphGroup);
		nodeRBGroup.addToParent(graphGroup);
		edgeRBGroup.addToParent(graphGroup);
		
		return graphGroup;
	}
	
	private static void addLightsToBranchGroup(BranchGroup bgroup) {
		Color3f nodeLightColor      = new Color3f(1.0f, 1.0f, 1.0f);
		BoundingSphere bounds       = new BoundingSphere(new Point3d(0,0,0), 100);
		Vector3f nodeLightDirection = new Vector3f(+1.0f, -1.0f, -1.0f);
		DirectionalLight nodeLight  = new DirectionalLight(nodeLightColor, nodeLightDirection);
		nodeLight.setInfluencingBounds(bounds);
		bgroup.addChild(nodeLight);
	}

	private static ReattachableBranchGroup createNodeGroup(CyNetworkView ourNetView, Point2D.Float[] srcBounds, double sizeOnScreen, HashMap rememberBranches) {
		// sourceBounds: the boundaries of the original source graph. We will want to scale the graph to lie in -1...+1 in XY
		CyNetwork   ourNetwork  = ourNetView.getNetwork();
		Appearance defaultAppear = DrawGraph3D.makeNodeAppearance();
		
		ReattachableBranchGroup nodeSphereGroup = new ReattachableBranchGroup();
		ReattachableBranchGroup nodeLabelGroup  = new ReattachableBranchGroup();	
		
		int     sphereDivisions = 4;
		float   sphereRadius    = 0.01f;
		int     primitiveFlags  = Primitive.ENABLE_GEOMETRY_PICKING + Primitive.GENERATE_NORMALS; // Primitive.GENERATE_TEXTURE_COORDS
		double  nodeZLoc        = 0.0;
		
		Point2D.Float translateAmt = getDrawingTranslation(srcBounds);
		double         scaleFactor = getDrawingScaleFactor(srcBounds, sizeOnScreen);
		
		Point3f rotateAroundPoint = new Point3f(0f, 0f, 0f);
		// Text labels will rotate around this point (local coords). (0,0,0) would rotate around the center of the text (this is a good default
				
		for (Iterator it = ourNetView.getNodeViewsIterator(); it.hasNext(); ) {
			NodeView nView = (NodeView)it.next();
			CyNode    node = (CyNode)nView.getNode();
			// Note that we flip the y-axis for the node's coordinates. This is 
			// because Cytoscape and the Java3D system have different coordinate systems.
			// When we move the graph around, translate it first, THEN apply the scale factor!
			Point3f nodePoint = new Point3f((float)(scaleFactor * (nView.getXPosition() + translateAmt.x)), 
										    (float)(-1*scaleFactor * (nView.getYPosition() + translateAmt.y)),
											(float) nodeZLoc);
			
			Transform3D  nodePositionT3D = new Transform3D();
			nodePositionT3D.setTranslation(new Vector3f(nodePoint));
			
			{ // NODE SPHERES
				TransformGroup tgroup = new TransformGroup();
				tgroup.setTransform(nodePositionT3D);
				Sphere nodeShape = new Sphere(sphereRadius, primitiveFlags, sphereDivisions, defaultAppear);
				//System.err.println("Added node at " + nodePoint + ".");
				tgroup.addChild(nodeShape);
				nodeSphereGroup.addChild(tgroup);
			}
			
			{ // LABELS FOR NODES
				String nodeCanonicalName = (String)(ourNetwork.getNodeAttributeValue(node, Semantics.CANONICAL_NAME));
				
				Text2D labelText2D = new Text2D(nodeCanonicalName, new Color3f(1.0f, 1.0f, 0.0f), "Helvetica", 36, 0);
				labelText2D.setRectangleScaleFactor(labelText2D.getRectangleScaleFactor() * 0.05f);
				//System.out.println("Scaling factor: " + texto.getRectangleScaleFactor());
								
				OrientedShape3D labelShape3D = new OrientedShape3D(labelText2D.getGeometry(), labelText2D.getAppearance(), OrientedShape3D.ROTATE_ABOUT_POINT, rotateAroundPoint);
				// labelShape3D: a "billboarded" label, always facing the viewer no matter what
				TransformGroup labelTGroup = new TransformGroup();
				labelTGroup.setTransform(nodePositionT3D); // same location as the node sphere...
				labelTGroup.addChild(labelShape3D);
				nodeLabelGroup.addChild(labelTGroup);
			}
			
		} // end for
		
		ReattachableBranchGroup bgroup = new ReattachableBranchGroup();
		nodeSphereGroup.addToParent(bgroup);
		nodeLabelGroup.addToParent(bgroup);
		if (rememberBranches != null) {
			rememberBranches.put(Topomap3D.kGRAPH_NODES, nodeSphereGroup);
			rememberBranches.put(Topomap3D.kGRAPH_NODE_LABELS, nodeLabelGroup);
		}
		return bgroup;
	} // end addNodesToBranchGroup
	
	private static double getDrawingScaleFactor(Point2D.Float[] srcBounds, double sizeOnScreen) {
		// Figure out how to scale the source graph so that it fits exactly onto the 3D drawing zone.
		double maxDistAcrossGraph = Math.max(srcBounds[ExtendedCyto.MAX_PT].x - srcBounds[ExtendedCyto.MIN_PT].x,
											 srcBounds[ExtendedCyto.MAX_PT].y - srcBounds[ExtendedCyto.MIN_PT].y);
		double scaleFactor = (sizeOnScreen / maxDistAcrossGraph);
		return scaleFactor;
	}
	
	private static Point2D.Float getDrawingTranslation(Point2D.Float[] srcBounds) {
		// Figure out how to translate the source graph (BEFORE SCALING) so that it fits exactly on the 3D drawing zone
		Point2D.Float midPt = new Point2D.Float((float) 0.5 * (srcBounds[ExtendedCyto.MIN_PT].x + srcBounds[ExtendedCyto.MAX_PT].x),
												(float) 0.5 * (srcBounds[ExtendedCyto.MIN_PT].y + srcBounds[ExtendedCyto.MAX_PT].y));
		return new Point2D.Float(-1 * midPt.x, -1 * midPt.y);
	}
	
	private static ReattachableBranchGroup createEdgeGroup(CyNetworkView ourNetView, Point2D.Float[] srcBounds, double sizeOnScreen) {
		CyNetwork theNetwork = ourNetView.getNetwork();
		
		Point2D.Float translateAmt = getDrawingTranslation(srcBounds);
		double         scaleFactor = getDrawingScaleFactor(srcBounds, sizeOnScreen);
		ReattachableBranchGroup bgroup = new ReattachableBranchGroup();
		
		final int numVertices = (2 * theNetwork.getEdgeCount());
		LineArray edgeGeom = new LineArray(numVertices, LineArray.COORDINATES); //|LineArray.COLOR_3
		
		double zLoc = 0.0; // <-- z-location where we draw the edges
		
		int   index = 0; // (the loop counter variable)
		for (Iterator it = theNetwork.edgesIterator(); it.hasNext(); ) {
			CyEdge edge = (CyEdge)it.next();
			CyNode sourceN = (CyNode)edge.getSource();
			CyNode targetN = (CyNode)edge.getTarget();
			
			NodeView sourceNV = ourNetView.getNodeView(sourceN);
			NodeView targetNV = ourNetView.getNodeView(targetN);
			
			// Below: Note that we flip the y-axis for the node's coordinates. This is 
			// because Cytoscape and the Java3D system have different coordinate systems.
			
			Point3d sourcePt = new Point3d(scaleFactor * (sourceNV.getXPosition() + translateAmt.x), 
										   -1*scaleFactor * (sourceNV.getYPosition() + translateAmt.y),
										   zLoc);
			edgeGeom.setCoordinate(index  , sourcePt);
			
			Point3d targetPt = new Point3d(scaleFactor * (targetNV.getXPosition() + translateAmt.x), 
										   -1*scaleFactor * (targetNV.getYPosition() + translateAmt.y),
										   zLoc);
			edgeGeom.setCoordinate(index+1, targetPt);
			
			index += 2;
			
			//System.err.println("Added edge from " + sourcePt + " to " + targetPt + ".");
			// A line array connects vertices like so:  0---1   2---3   4---5   6---7
		} // end for
		TransformGroup tgroup = new TransformGroup();
		Transform3D       t3d = new Transform3D();
		t3d.setTranslation(new Vector3d(0, 0, 0));
		tgroup.setTransform(t3d);
		
		Shape3D edgeShape = new Shape3D(edgeGeom, makeLineAppearance());
		//edgeShape.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE); // Note that only one capability bit may be set per invocation
		tgroup.addChild(edgeShape);
		
		bgroup.addChild(tgroup);
		return bgroup;
	} // end addEdgesToBranchGroup
	
	
	private static Appearance makeLineAppearance() {
		Appearance edgeAppearance = new Appearance();
		LineAttributes lineAttrib = new LineAttributes();
		lineAttrib.setLineWidth(1.0f);
		lineAttrib.setLinePattern(LineAttributes.PATTERN_SOLID);//0xFFFF;
			lineAttrib.setLineAntialiasingEnable(true);
			edgeAppearance.setLineAttributes(lineAttrib);
			
			TransparencyAttributes transparAttrib = new TransparencyAttributes();
			transparAttrib.setTransparency(0.5f); // 0.0 = solid, 1.0 = completely transparent
			transparAttrib.setTransparencyMode(TransparencyAttributes.BLENDED);
		//transparAttrib.setTransparencyMode(TransparencyAttributes.NICEST);
		//transparAttrib.setTransparencyMode(TransparencyAttributes.SCREEN_DOOR);
		//transparAttrib.setTransparencyMode(TransparencyAttributes.NONE);
			edgeAppearance.setTransparencyAttributes(transparAttrib);
			
			return edgeAppearance;
	}
	
	private static Appearance makeNodeAppearance() {
		Appearance nodeAppear = new Appearance();
		PolygonAttributes polyAttrib = new PolygonAttributes();
		polyAttrib.setCullFace(PolygonAttributes.CULL_BACK); // cull back-facing polygons
		polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_FILL); // cull back-facing polygons
		nodeAppear.setPolygonAttributes(polyAttrib);
		nodeAppear.setMaterial(new Material(nodeDefaultColor, nodeDefaultColor, nodeDefaultColor, nodeDefaultColor, 1.0f));
		// material contructor: Material(Color3f ambientColor, Color3f emissiveColor, Color3f diffuseColor, Color3f specularColor, float shininess) 
		return nodeAppear;
	}
	
}
