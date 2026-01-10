//  Topomap3D.java
//  Topomap3D
//
//  Created by Alex Williams on 9/21/05.


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

import java.util.HashMap;

import java.applet.Applet;
import java.awt.BorderLayout;
import java.awt.Frame;
import java.awt.event.*;
import java.awt.GraphicsConfiguration;
import com.sun.j3d.utils.applet.MainFrame; 
import com.sun.j3d.utils.universe.*;
import com.sun.j3d.utils.geometry.Text2D;

import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.behaviors.keyboard.*;
import com.sun.j3d.utils.behaviors.mouse.*;
import com.sun.j3d.utils.behaviors.sensor.SensorGnomonEcho;

import java.util.Enumeration;
import java.util.Locale;
import java.util.ResourceBundle;
import java.util.BitSet;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

import javax.swing.*;

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

import org.newdawn.common.hud.*;

//import com.apple.eawt.*;

// http://download.java.net/media/java3d/javadoc/1.4.0-latest/com/sun/j3d/utils/behaviors/mouse/MouseBehavior.html

public class Topomap3D extends JPanel {
	
	public static final String kTERRAIN = "Terrain";
	public static final String kGRAPH_ENTIRE = "Graph (Entire)";
	public static final String kGRAPH_NODES = "Graph (Nodes)";
	public static final String kGRAPH_EDGES = "Graph (Edges)";
	public static final String kGRAPH_NODE_LABELS = "Graph Node Labels";
	
	private static boolean STANDALONE_MODE = false; // is it running outside of Cytoscape ("standalone mode") or inside it?
	
	private static final double kBACK_CLIPPING_DISTANCE  = 10.0; // how far away things are no longer drawn
	private static final double kFRONT_CLIPPING_DISTANCE = 0.005; // how close until things are no long drawn
	private static final double kBOUNDING_SPHERE_RADIUS = 100; // needs to be large enough to include everything in the scene
	private static final long   kMILLISECONDS_BETWEEN_UPDATES = 25; // smaller values = more updates
	private static final int    kINITIAL_WINDOW_SIZE = 800; // in pixels--only works for the standalone mode, though
	
	
	private float   topomapVerticalScale;
	public float getTopomapVerticalScale() { return this.topomapVerticalScale; }
	public void setTopomapVerticalScale(float x) { this.topomapVerticalScale = x; }
	
	private boolean showTopomap;
	public boolean getShowTopomap() { return this.showTopomap; }
	public void setShowTopomap(boolean x) { this.showTopomap = x; }
	
	private boolean showGraphNodes, showGraphEdges, showGraphLabels;
	public boolean getShowGraphNodes() {  return this.showGraphNodes; }
	public boolean getShowGraphEdges() {  return this.showGraphEdges; }
	public boolean getShowGraphLabels() { return this.showGraphLabels; }
	public void setShowGraphNodes(boolean x) {  this.showGraphNodes = x; }
	public void setShowGraphEdges(boolean x) {  this.showGraphEdges = x; }
	public void setShowGraphLabels(boolean x) { this.showGraphLabels = x; }
	
	private HUD hud;
	private HUDLegend hudLegend;

	private CustomCamera cam;
	public CustomCamera getCamera() { return cam; }
	
	//private Transform3D cameraPoint, cameraLook; // location of camera and direction in which it's looking
	
	private Font font = new Font("serif", Font.ITALIC+Font.BOLD, 36);
	
	private static final JMenuBar mainMenuBar = new JMenuBar();
	
	private Canvas3D canvas;
	private Shape3D  terrain;
	public Shape3D getTerrainShape3D() { return terrain; }
	
	public BranchGroup graphBranch; // the ENTIRE 2D graph (nodes & edges)
	
	public HashMap branches = new HashMap(); // a hashmap that keeps track of all the ReattachableBranchGroups
	
	private SimpleUniverse univ;
	public SimpleUniverse getUniverse() { return univ; }
	
	private CustomHandleKeys keyBehavior;
	private CustomHandleMouse mouseBehavior;
	
	private float fieldOfView = 120.0f;
	
	// No main function anymore!
	
	Topomap3D(Density density) {
		 super();
		 //this.getContentPane().setLayout(new BorderLayout());
		 GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
		 canvas = new Canvas3D(config);
		 
		 this.topomapVerticalScale = 1.0f;
		 this.showTopomap = true;
		 this.showGraphNodes = this.showGraphEdges = this.showGraphLabels = true;
		 
		 this.setLayout(new BorderLayout());
		 this.add(canvas, BorderLayout.CENTER);
			
		 this.univ = new SimpleUniverse(canvas);
		 //this.univ.getViewer().getView().setTransparencySortingPolicy(View.TRANSPARENCY_SORT_GEOMETRY);
		 
		 BranchGroup scene = createSceneGraph(density);
		 
		 this.univ.addBranchGraph(scene);
		 
		 this.univ.getViewer().getView().setFieldOfView(this.fieldOfView);
		 this.univ.getViewer().getView().setBackClipDistance(kBACK_CLIPPING_DISTANCE); // how far away things are no longer drawn
		 this.univ.getViewer().getView().setFrontClipDistance(kFRONT_CLIPPING_DISTANCE); // how close until things are no long drawn
		 this.univ.getViewer().getView().setSceneAntialiasingEnable(true);
		 
		 //setSize(900, 900);
		 
		 TransformGroup vpTG = this.univ.getViewingPlatform().getViewPlatformTransform();
		 
		 this.cam = new CustomCamera();
		 //this.univ.getViewingPlatform().setNominalViewingTransform();
		 this.cam.moveCameraTo(new Point3d(0, 0, 5),
							   new Point3d(0, 0, 0));
		 this.redoView();    // used for initial position 
		 
		 this.setVisible(true);
	}
	
	public void redoView() {
		// Viewpoint...
		(this.univ.getViewingPlatform().getViewPlatformTransform()).setTransform(this.cam.getTotalTransformation());
	}
		
	public void redoBranches() {
		// Redetermines which branches are visible and which ones are not, based on the booleans in this function.
		try {
			((ReattachableBranchGroup)branches.get(kTERRAIN)).setAttach(this.showTopomap);
			((ReattachableBranchGroup)branches.get(kGRAPH_NODES)).setAttach(this.showGraphNodes);
			((ReattachableBranchGroup)branches.get(kGRAPH_EDGES)).setAttach(this.showGraphEdges);
			((ReattachableBranchGroup)branches.get(kGRAPH_NODE_LABELS)).setAttach(this.showGraphLabels);
			
		} catch (Exception e) {
			// Report any problems, but don't exit the program.
			e.printStackTrace();
			System.err.println(e);
		}
		
		this.redrawHUD();
		this.redoView(); // And make sure to redraw things, too
	}
	
	public void redrawHUD() {
		hudLegend.draw();
		//HUDComponent[] components = this.hud.getHUDComponents();
		//for (int i = 0; i < components.length; i++) {
		//	((HUDLegend)components[i]).paint(components[i]);
		//}
	}
	
	public void testMark(Topomap3D obj, float r, float g, float b, int num) {
		LandscapeShape3D lg = (LandscapeShape3D)(obj.getTerrainShape3D());
		lg.randomizeTerrain(); // randomize terrain
		lg.makeLandscapeGeometry();
		System.out.println("Remaking terrain...");
		
		//lg.setGeometry(lg.makeLandscapeGeometry(null));
		//TriangleStripArray t2 = (TriangleStripArray) obj.getTerrainShape3D().getGeometry();
		//t2.setCoordinate(num, new Point3f(0,0, 1));
		//t2.setColor(num, new Color3f(r, g, b));
	}
	
	public BranchGroup createSceneGraph(Density density) {
		
		BranchGroup contentRoot = new BranchGroup();
		contentRoot.setCapability(BranchGroup.ALLOW_DETACH);
		contentRoot.setCapability(Group.ALLOW_CHILDREN_READ);
		contentRoot.setCapability(Group.ALLOW_CHILDREN_WRITE);
		contentRoot.setCapability(Group.ALLOW_CHILDREN_EXTEND);
		
		// Heirarchy of geometry / transform groups:
		//    contentRoot
		//    |         \
		//    |          background 
		// interactiveTG
		//    |        \  
		// terrainTG   graphTG
		//    |        (this.graphBranch)
		//    |
		// (this.terrain)
		
		// everything you can interact with in a normal way goes into this TG		
		TransformGroup interactiveTG = new TransformGroup();
		interactiveTG.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		interactiveTG.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
		interactiveTG.setCapability(Group.ALLOW_CHILDREN_READ);
		interactiveTG.setCapability(Group.ALLOW_CHILDREN_WRITE);
		interactiveTG.setCapability(Group.ALLOW_CHILDREN_EXTEND);
		
		// TERRAIN
		Transform3D terrainTransform = new Transform3D();
		// Put the terrain into a cube from (-0.5, -0.5, 0) to (+0.5, +0.5, +1)
		terrainTransform.setTranslation(new Vector3f(-0.5f, -0.5f, 0));		
		TransformGroup terrainTG = new TransformGroup(terrainTransform); 
		terrainTG.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		terrainTG.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
		this.terrain = new LandscapeShape3D(density);
		terrainTG.addChild( (Shape3D)this.terrain );
		ReattachableBranchGroup terrainRBG = new ReattachableBranchGroup();
		branches.put(kTERRAIN, terrainRBG); // add the terrain RBG to the list of all the reattachable branches
		terrainRBG.addChild(terrainTG);
		terrainRBG.addToParent(interactiveTG); // adds terrainRBG as a child of interactiveTG
		
		// GRAPH
		
		if (!STANDALONE_MODE) {
			// Only call Cytoscape methods if we are NOT running stand-alone
			TransformGroup graphTG = new TransformGroup();
			graphTG.setCapability(Group.ALLOW_CHILDREN_READ);
			graphTG.setCapability(Group.ALLOW_CHILDREN_WRITE);
			graphTG.setCapability(Group.ALLOW_CHILDREN_EXTEND);
			
			this.graphBranch = DrawGraph3D.newGraph(Cytoscape.getCurrentNetworkView(), density.getNonBorderFraction(), branches);
			//this.graphBranch.setCapability(BranchGroup.ALLOW_DETACH);
			this.graphBranch.setCapability(Group.ALLOW_CHILDREN_READ);
			this.graphBranch.setCapability(Group.ALLOW_CHILDREN_WRITE);
			this.graphBranch.setCapability(Group.ALLOW_CHILDREN_EXTEND);
			
			graphTG.addChild(this.graphBranch);
			interactiveTG.addChild(graphTG);   // Add the nodes & edges
		}
		
		BoundingSphere bounds = new BoundingSphere(new Point3d(), kBOUNDING_SPHERE_RADIUS);
		// Effectively infinite boundaries, since our geometry usually goes from -1 to 1
		
		// Our custom keyboard input handler
		keyBehavior = new CustomHandleKeys(interactiveTG, this);
		keyBehavior.setSchedulingBounds(bounds);
		interactiveTG.addChild(keyBehavior);
		
		// Our custom mouse input handler
		//mouseBehavior = new CustomHandleMouse(interactiveTG, this);
		//mouseBehavior.setSchedulingBounds(bounds);
		//contentRoot.addChild(mouseBehavior);
		
		MouseRotate myMouseRotate = new MouseRotate(); //MouseBehavior.INVERT_INPUT); 
		myMouseRotate.setTransformGroup(interactiveTG); 
		myMouseRotate.setSchedulingBounds(bounds); 
		interactiveTG.addChild(myMouseRotate);

		CustomTimer timer = new CustomTimer(kMILLISECONDS_BETWEEN_UPDATES, this);
		timer.setSchedulingBounds(bounds);
		interactiveTG.addChild(timer);

		contentRoot.addChild(interactiveTG);
		
		Background background = new Background(0.0f, 0.0f, 0.0f);
		background.setApplicationBounds(bounds);
		contentRoot.addChild(background);
		
		canvas.getView().setFrontClipPolicy(View.VIRTUAL_EYE);
		canvas.getView().setTransparencySortingPolicy(View.TRANSPARENCY_SORT_GEOMETRY);
		
		// HUD things
		int numLayers = 2;
		float baseDepth = 0.1f;
		this.hud = new PhysicalHUD(numLayers, this.univ, baseDepth);
		this.hudLegend = new HUDLegend(this.hud, this);
		//TestDisplayArea testDisplayArea = new TestDisplayArea(hud);
		//TestInputArea testInputArea = new TestInputArea(hud, testDisplayArea);
		//TestGrid testGrid = new TestGrid(hud);
		//TestLabel tl   = new TestLabel(hud);
		
		contentRoot.compile(); // optimize the scene graph (optional step)
		
		return contentRoot;
	}
	
	public void handleTimeUpdate() {
		keyBehavior.handleTimeUpdate();
		//System.err.println("Timer called!\n");
		this.redoView();
	}
	
	public void quit() {	System.exit(0);		}
	
	public static void main(String args[]) {
		// This is called ONLY when the program is run as a stand-alone application.
		// It is NOT EVER called when we run the program from inside Cytoscape.
		STANDALONE_MODE = true; // we are running stand-alone
		JFrame appFrame = new JFrame("Title Bar Text");
		JPanel mainPanel = new Topomap3D(null); //new TopomapPanel(null);
		appFrame.setSize(kINITIAL_WINDOW_SIZE, kINITIAL_WINDOW_SIZE);
		appFrame.setDefaultCloseOperation(appFrame.EXIT_ON_CLOSE);
		appFrame.getContentPane().add(mainPanel);
		appFrame.show();
		System.err.println("\nRunning the Topomap program in standalone mode.\n");
	}
	
}


// Extra / old code

//SensorGnomonEcho unitScaleAxes = new SensorGnomonEcho(new Transform3D(), 0.05, 0.5, false); // one unit "diameter" from the tips of the arms

//		MouseZoom myMouseZoom = new MouseZoom();
//		myMouseZoom.setTransformGroup(interactiveTG); 
//		myMouseZoom.setSchedulingBounds(bounds); 
//		contentRoot.addChild(myMouseZoom); 

//		MouseTranslate myMouseTranslate = new MouseTranslate();
//		myMouseTranslate.setTransformGroup(interactiveTG); 
//		myMouseTranslate.setSchedulingBounds(bounds); 
//		contentRoot.addChild(myMouseTranslate);
