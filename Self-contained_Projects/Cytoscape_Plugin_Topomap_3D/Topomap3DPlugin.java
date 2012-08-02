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

import java.util.*;
import java.awt.event.ActionEvent;
import javax.swing.JOptionPane;

import giny.model.Node;
import giny.view.NodeView;

import cytoscape.plugin.CytoscapePlugin;
import cytoscape.util.CytoscapeAction;
import cytoscape.Cytoscape;
import cytoscape.CyNetwork;
import cytoscape.CyNode;
import cytoscape.CyEdge;
import cytoscape.view.CyNetworkView;
import cytoscape.data.Semantics;
import cytoscape.data.readers.GraphReader;

import java.io.FileWriter;

import java.util.ArrayList;
import java.awt.geom.Point2D;

// temporary imports, only for testing
import javax.swing.*;
import java.awt.*;

public class Topomap3DPlugin extends CytoscapePlugin {
	
	public String describe() {
		return "Test plugin to try Java 3D";
		
	}
	
	public Topomap3DPlugin() {		
		// Add a menu option for running this plugin in the main cytoscape view window
		//TopoAction la = new TopoAction();
		//la.setPreferredMenu("Plugins");
		//Cytoscape.getDesktop().getCyMenus().addAction(la);// Add fullscreen mode menu option
		
		CytoscapeAction ela = new CytoscapeAction("Topomap") {
			public void actionPerformed(ActionEvent ae) {
				String title = "Topomap 3D Graph";
				int numBorderSpaces = 6;
				int numTotalSpaces   = 180;
				Density d = new Density(getNodePositionsFromCurrentView(),
										numBorderSpaces, numTotalSpaces, 
										Density.makeQuadraticKernel(5, -0.04f, 0, 1));
				CyNetworkView view = new TopoNetworkView(Cytoscape.getCurrentNetwork(), title, d);
				Topomap3DPlugin.addCustomNetworkView(view);
			}
		};
		ela.setAcceleratorCombo(java.awt.event.KeyEvent.VK_T, java.awt.event.ActionEvent.CTRL_MASK);
		ela.setPreferredMenu("Plugins");
		Cytoscape.getDesktop().getCyMenus().addAction(ela); // Add windowed mode menu option
		
		BoundingBoxAction ba = new BoundingBoxAction();
		ba.setPreferredMenu("Plugins");
		Cytoscape.getDesktop().getCyMenus().addAction(ba); // Add another menu option
		
		MultiLevelLayoutAction mla = new MultiLevelLayoutAction("Multi-level Layout");
		mla.setPreferredMenu("Plugins");
		mla.setAcceleratorCombo(java.awt.event.KeyEvent.VK_M, java.awt.event.ActionEvent.CTRL_MASK);
		Cytoscape.getDesktop().getCyMenus().addAction(mla); // Add another menu option
		
		CytoscapeAction reclusterAction = new CytoscapeAction("Multi-level Layout: Re-cluster around centroids") {
			public void actionPerformed(ActionEvent ae) {
				System.err.println("Starting recentering...");
				MultiLevelLayoutAction.recenterAroundCurrentCentroids(Cytoscape.getCurrentNetworkView());
				System.err.println("Done recentering.");
			}
		};
		reclusterAction.setPreferredMenu("Plugins");
		reclusterAction.setAcceleratorCombo(java.awt.event.KeyEvent.VK_K, java.awt.event.ActionEvent.CTRL_MASK);
		Cytoscape.getDesktop().getCyMenus().addAction(reclusterAction); // Add another menu option
		
		DensityAction da = new DensityAction();
		da.setPreferredMenu("Plugins");
		Cytoscape.getDesktop().getCyMenus().addAction(da); // Add another menu option
		
		/*
		CytoscapeAction ea = new CytoscapeAction("Throw Exception") {
			public void actionPerformed(ActionEvent ae) {
				Object o = null;
				o.toString();
			}
		};
		ea.setPreferredMenu("Plugins");
		Cytoscape.getDesktop().getCyMenus().addAction(ea);  // Add another menu option
		*/
	}

	// this assumes that the view is for a network that has already been
	// loaded, which is far as i can tell will always be true.
	static public CyNetworkView addCustomNetworkView(CyNetworkView view) {
		CyNetwork net = view.getNetwork();
		
		// at this point in time, depend on an existing 2D view
		if (!Cytoscape.viewExists(net.getIdentifier())) {return null;}
		
		System.out.println("attempting to destroy previous view");
		Cytoscape.destroyNetworkView(net);
		
		// this is modeled after Cytoscape.createNetworkView() without
		// forcing the view to be a PhoebeNetworkView
		view.setIdentifier(net.getIdentifier());
		Cytoscape.getNetworkViewMap().put(net.getIdentifier(), view);
		view.setTitle(net.getTitle());
		
		if (net.getClientData("GML") != null) {
			((GraphReader) net.getClientData("GML")).layout(view);
		}
		
		Cytoscape.firePropertyChange(cytoscape.view.CytoscapeDesktop.NETWORK_VIEW_CREATED, null, view);
		
		view.redrawGraph(false, false);
		
		return view;
	}
	
	// Generate an array of the point locations of all the nodes in a view
	// Warning: I'm not sure if the view is allowed to change while
	// plugins are running... if so, then the weird iterator/array
	// conversion thing will go awry.
	static public Point2D[] getNodePositionsFromView(CyNetworkView view) {
		if (view==null) {return null;}
		
		Point2D[] result = new Point2D[view.getNodeViewCount()];
		
		Iterator iter = view.getNodeViewsIterator();
		for (int i = 0; i < result.length; ++i) {
			NodeView nv = (NodeView)iter.next();
			result[i] = new Point2D.Double(nv.getXPosition(), nv.getYPosition());
		}
		return result;
	}
	
	static public Point2D[] getNodePositionsFromCurrentView() {
		return getNodePositionsFromView(Cytoscape.getCurrentNetworkView());
	}
	
	public class TopoAction extends CytoscapeAction {
		// Define a menu option for running this plugin in a separate window
		public TopoAction() {super("Topomap (Separate window)");}
        
		public void actionPerformed(ActionEvent ae) {
			Density d = new Density(getNodePositionsFromCurrentView(),
									9, 150,
									Density.makeQuadraticKernel(5,-0.04f,0,1));
			
			new Topomap3D(d);
			String message = "Running Topomap in a separate window...";
			System.out.println(message);
			//JOptionPane.showMessageDialog(Cytoscape.getDesktop(), message);
		}
	}
	
	public class BoundingBoxAction extends CytoscapeAction {
		public BoundingBoxAction() {super("Bounding Box");}
		public void actionPerformed(ActionEvent ae) {
			Point2D[] b = Density.getBoundsOfPoints(getNodePositionsFromCurrentView());
			
			if (b == null) {return;}
			
			// Output bounding box
			try {
				FileWriter out = new FileWriter("boundingbox.tab", false);
				out.write("xrange\t" + b[0].getX() + "\t" + b[1].getX() + "\n");
				out.write("yrange\t" + b[0].getY() + "\t" + b[1].getY() + "\n");
				out.close();
			} catch (java.io.IOException ioe) {
				System.err.println("Couldn't open bounding box output file");
			}
		} // void actionPerformed()
	} // class densityAction
	
	
	public class DensityAction extends CytoscapeAction {
		public DensityAction() {super("Write Density File");} 
		
		public void actionPerformed(ActionEvent ae) {
			Density d = new Density(getNodePositionsFromCurrentView());
			
			try {
				FileWriter out = new FileWriter("density.tab", false);
				out.write(d.toString());
				out.close();
			} catch (java.io.IOException ioe) {
				System.err.println("Couldn't write file 'density.tab'");
			}
		}
	}
} // class Topomap3DPlugin
