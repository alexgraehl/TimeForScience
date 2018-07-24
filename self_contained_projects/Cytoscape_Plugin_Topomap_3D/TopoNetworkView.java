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


import java.util.Map;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Iterator;
import java.util.Arrays;

import java.awt.Component;
import java.awt.Paint;

import cytoscape.*;
import cytoscape.layout.*;
import cytoscape.view.*;
import cytoscape.visual.*;

import giny.view.*;
import giny.model.*;

import coltginy.ColtGraphViewModel;

public class TopoNetworkView extends ColtGraphViewModel implements CyNetworkView {
  protected CyNetwork network;
  protected String title;
  protected String identifier;

  protected boolean isVizMapping = false;

  protected Map clientData;

  protected Component viewComponent;

  // selection event handler
  protected boolean nodeSelection;
  protected boolean edgeSelection;

  public TopoNetworkView(CyNetwork network, String title, Density d) {
    super(title, network);
    this.network = network;
    this.title = title;

    this.clientData = new HashMap();
    
    System.out.println("Creating new TopoNetworkView");
    
    this.viewComponent = new Topomap3D(d); // new TopomapPanel(d);
  }

  public CyNetwork getNetwork () { return network; }
  public void setTitle ( String title ) { this.title = title; }
  public String getTitle () { return title; }
  
  public void redrawGraph( boolean layout, boolean vizmap ) {
    
  }

  public CyNetworkView getView () { return (CyNetworkView)this; }
 
  public VisualMappingManager getVizMapManager() {
    return null;
  }

  public cytoscape.visual.ui.VizMapUI getVizMapUI() {
    return null; // this is the same as PhoebeNetworkViewer
  }
  
  public void toggleVisualMapperEnabled() {
    // Don't do anything, yet
  }

  public void setVisualMapperEnabled ( boolean state ) {
    // Don't do anything, yet
  }

  public boolean getVisualMapperEnabled () { return isVizMapping; }

  //--------------------//
  // Network Client Data
  
  public void putClientData ( String data_name, Object data ) {
    clientData.put(data_name, data);
  }

  public Collection getClientDataNames () {
    return clientData.keySet();
  }    

  public Object getClientData ( String data_name ) {
    return clientData.get(data_name);
  }
    

  public boolean setSelected ( CyNode[] nodes ) {
    return setSelected(convertToViews(nodes));
  }
  public boolean setSelected ( NodeView[] node_views ) {
    for (int i = 0; i < node_views.length; ++i) {
      node_views[i].select();
    }
    return true;
  }

  public boolean setSelected ( CyEdge[] edges ) {
    return setSelected(convertToViews(edges));
  }
  public boolean setSelected ( EdgeView[] edge_views ) {
    for (int i = 0; i < edge_views.length; ++i) {
      edge_views[i].select();
    }
    return true;
  }
 

  /**
   * Applies the given edge to the default vizmapper
   */
  public boolean applyVizMap ( CyEdge edge ) {
    return applyVizMap(getEdgeView(edge));
  }

  /**
   * Applies the given edge to the default vizmapper
   */
  public boolean applyVizMap (EdgeView edge_view) {
    return
      applyVizMap(edge_view, 
                  (VisualStyle)getClientData(CytoscapeDesktop.VISUAL_STYLE));
  }

  /**
   * Applies the given edge to the given vizmapper
   */
  public boolean applyVizMap ( CyEdge edge, VisualStyle style ) {
    return applyVizMap(getEdgeView(edge), style);
  }

  /**
   * Applies the given edge to the given vizmapper
   */
  public boolean applyVizMap ( EdgeView edge_view, VisualStyle style ) {
    VisualStyle old_style = Cytoscape.getDesktop().setVisualStyle( style );
    Cytoscape.getDesktop().getVizMapManager().vizmapEdge( edge_view, this );
    Cytoscape.getDesktop().setVisualStyle( old_style );
    return true;
  }

  /**
   * Applies the given node to the default vizmapper
   */
  public boolean applyVizMap ( CyNode node ) {
    return
      applyVizMap(node, 
                  (VisualStyle)getClientData(CytoscapeDesktop.VISUAL_STYLE));
  }

  /**
   * Applies the given node to the default vizmapper
   */
  public boolean applyVizMap ( NodeView node_view ) {
    return 
      applyVizMap(node_view, 
                  (VisualStyle)getClientData(CytoscapeDesktop.VISUAL_STYLE));
  }

  /**
   * Applies the given node to the given vizmapper
   */
  public boolean applyVizMap ( CyNode node, VisualStyle style ) {
    return applyVizMap(getNodeView(node), style);
  }

  /**
   * Applies the given node to the given vizmapper
   */
  public boolean applyVizMap ( NodeView node_view, VisualStyle style ) {
    VisualStyle old_style = Cytoscape.getDesktop().setVisualStyle(style);
    Cytoscape.getDesktop().getVizMapManager().vizmapNode(node_view, this);
    Cytoscape.getDesktop().setVisualStyle(old_style);
    return true;
  }

  public void applyVizmapper ( VisualStyle style ) {
    VisualStyle old_style = Cytoscape.getDesktop().setVisualStyle( style );
    redrawGraph( false, true );
    Cytoscape.getDesktop().setVisualStyle(old_style);
  }
    
  /**
   * Applies the given layout to the entire CyNetworkView
   */
  public void applyLayout ( LayoutAlgorithm layout ) {
    layout.doLayout();
  }

  /**
   * Applies the given layout to the entire CyNetworkView,
   * but locks the given Nodes and Edges in place
   */
  public void applyLockedLayout ( LayoutAlgorithm layout, CyNode[] nodes, CyEdge[] edges ) {
    layout.lockNodes( convertToViews( nodes ) );
    layout.doLayout();
  }

  /**
   * Applies the  given layout to only the given Nodes and Edges
   */
  public void applyLayout ( LayoutAlgorithm layout, CyNode[] nodes, CyEdge[] edges) {
    layout.lockNodes( getInverseViews( convertToViews( nodes ) ) );
    layout.doLayout();
  }


  /**
   * Applies the given layout to the entire CyNetworkView,
   * but locks the given NodeViews and EdgeViews in place
   */
  public void applyLockedLayout ( LayoutAlgorithm layout, CyNodeView[] nodes, CyEdgeView[] edges ) {
    layout.lockNodes(  nodes );
    layout.doLayout();
  }

  /**
   * Applies the  given layout to only the given NodeViews and EdgeViews
   */
  public void applyLayout ( LayoutAlgorithm layout, CyNodeView[] nodes, CyEdgeView[] edges ) {
    layout.lockNodes( getInverseViews( nodes ) );
    layout.doLayout();
  }

  /**
   * Applies the given layout to the entire CyNetworkView,
   * but locks the given Nodes and Edges in place
   */
  public void applyLockedLayout ( LayoutAlgorithm layout, int[] nodes, int[] edges ) {
    layout.lockNodes( convertToNodeViews( nodes ) );
    layout.doLayout();
  }

  /**
   * Applies the  given layout to only the given Nodes and Edges
   */
  public void applyLayout ( LayoutAlgorithm layout, int[] nodes, int[] edges ) {

    layout.lockNodes( getInverseViews( convertToNodeViews( nodes ) ) );
    layout.doLayout();
  }

  //**********************
  // Helper functions
  //**********************
  protected NodeView[] getInverseViews ( NodeView[] given ) {
    NodeView[] inverse = new NodeView[ getNodeViewCount() - given.length ];
    List node_views = getNodeViewsList();
    int count = 0;
    Iterator i = node_views.iterator();
    Arrays.sort( given );
    while ( i.hasNext() ) {
      NodeView view = ( NodeView )i.next();
      if ( Arrays.binarySearch( given, view ) < 0 ) {
        // not a given, add
        inverse[count] = view;
        count++;
      }
    }
    return inverse;
  }

  protected EdgeView[] getInverseViews ( EdgeView[] given ) {
    EdgeView[] inverse = new EdgeView[ getEdgeViewCount() - given.length ];
    List edge_views = getEdgeViewsList();
    int count = 0;
    Iterator i = edge_views.iterator();
    Arrays.sort( given );
    while ( i.hasNext() ) {
      EdgeView view = ( EdgeView )i.next();
      if ( Arrays.binarySearch( given, view ) < 0 ) {
        // not a given, add
        inverse[count] = view;
        count++;
      }
    }
    return inverse;
  }

  protected NodeView[] convertToViews (CyNode[] nodes) {
    NodeView[] views = new NodeView[nodes.length];
    for (int i = 0; i < nodes.length; ++i) {
      views[i] = getNodeView(nodes[i]);
    }
    return views;    
  }

  protected EdgeView[] convertToViews (CyEdge[] edges) {
    EdgeView[] views = new EdgeView[edges.length];
    for (int i = 0; i < edges.length; ++i) {
      views[i] = getEdgeView(edges[i]);
    }
    return views;    
  }


  protected NodeView[] convertToNodeViews (int[] nodes) {
    NodeView[] views = new NodeView[nodes.length];
    for (int i = 0; i < nodes.length; ++i) {
      views[i] = getNodeView(nodes[i]);
    }
    return views;    
  }

  protected EdgeView[] convertToEdgeViews (int[] edges) {
    EdgeView[] views = new EdgeView[edges.length];
    for (int i = 0; i < edges.length; ++i) {
      views[i] = getEdgeView(edges[i]);
    }
    return views;    
  }


  //************************************************************
  //************************************************************
  //
  // giny.view.GraphView methods
  //
  //************************************************************
  //************************************************************

  public GraphPerspective getGraphPerspective() {return network;}

  //----------------------------------------//
  // Event Handlers
  //----------------------------------------//
  public boolean nodeSelectionEnabled() {return nodeSelection;}
  public boolean edgeSelectionEnabled() {return edgeSelection;}
  public void enableNodeSelection () {nodeSelection = true;}
  public void disableNodeSelection () {nodeSelection = false;}
  public void enableEdgeSelection () {edgeSelection = true;}
  public void disableEdgeSelection () {edgeSelection=false;}


  public int[] getSelectedNodeIndices() {return null;}
  public List getSelectedNodes() {return null;}

  public int[] getSelectedEdgeIndices() {return null;}
  public List getSelectedEdges() {return null;}

  public void addGraphViewChangeListener(GraphViewChangeListener l) {}
  public void removeGraphViewChangeListener(GraphViewChangeListener l) {}

  public void setBackgroundPaint(Paint paint) {}
  public Paint getBackgroundPaint() {return null;}

  public Component getComponent() {
    return viewComponent;
  }
  
  public NodeView addNodeView(int node_index) {
    return null;
    // !!! implement before running?
  }
  public NodeView addNodeView(String class_name, int node_index) {
    return null;
    // !!! implement before running?
  }
  public NodeView addNodeView ( int node_index, NodeView replacement ) {
    return null;
  }

  public EdgeView addEdgeView(int edge_index) {
    return null;
    // !!! implement before running?
  }

  public EdgeView addEdgeView(String class_name, int edge_index) {
    return null;
    // !!! implement before running?
  }


  public NodeView removeNodeView(NodeView nodeview) {
    return null;
    // !!! implement before running?
  }
  public NodeView removeNodeView(Node node) {
    return removeNodeView(getNodeView(node));
  }
  public NodeView removeNodeView(int node) {
    return removeNodeView(getNodeView(node));
  }
    

  public EdgeView removeEdgeView(EdgeView edgeview) {
    return null;
    // !!! implement before running?
  }
  public EdgeView removeEdgeView(Edge edge) {
    return removeEdgeView(getEdgeView(edge));
  }
  public EdgeView removeEdgeView(int edge) {
    return removeEdgeView(getEdgeView(edge));
  }

  public String getIdentifier() {return identifier;}
  public void setIdentifier(String identifier) {this.identifier = identifier;}

  public double getZoom() { return 1.0; }
  public void setZoom(double zoom) {}

  public void fitContent() {}

  public void updateView() {} // this isn't implemented in Phoebe either

  public RootGraph getRootGraph() { 
    return getGraphPerspective().getRootGraph();
  }
  
  public List getNodeViewsList() {return null;}
  public Iterator getNodeViewsIterator() {return null;}
  public int getNodeViewCount() {return 0;}
  public NodeView getNodeView(Node node) {return null;}
  public NodeView getNodeView(int node) {return null;}
  public int nodeCount() {return 0;}
  
  public List getEdgeViewsList() {return null;}
  public List getEdgeViewsList(Node a, Node b) {return null;}
  public List getEdgeViewsList(int a, int b, boolean include_undirected) {
    return null;
  }
  public Iterator getEdgeViewsIterator() {return null;}
  public int getEdgeViewCount() {return 0;}
  public EdgeView getEdgeView(Edge edge) {return null;}
  public EdgeView getEdgeView(int edge) {return null;}
  public int edgeCount() {return 0;}

  
  public boolean hideGraphObject(Object object) {return false;}
  public boolean showGraphObject(Object object) {return false;}
  public boolean showGraphObjects(List objects) {return false;}
  public boolean hideGraphObjects(List objects) {return false;}

  public Object[] getContextMethods(String class_name,
                                    boolean plus_superclass) {return null;}
  public Object[] getContextMethods(String class_name,
                                    Object[] methods) {return null;}
  public boolean addContextMethod(String class_name,
                                  String method_class_name,
                                  String method_name,
                                  Object[] args,
                                  ClassLoader loader) {return false;}

}