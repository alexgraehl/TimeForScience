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

import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;
import java.awt.GraphicsConfiguration;
import com.sun.j3d.utils.applet.MainFrame; 
import com.sun.j3d.utils.universe.*;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.behaviors.keyboard.*;
import com.sun.j3d.utils.behaviors.mouse.*;

import java.util.Enumeration;
//import java.util.Locale;
//import java.util.ResourceBundle;
import java.util.BitSet;


public class CustomHandleMouse extends Behavior {

	private BitSet mouseDownMap = new BitSet(); // a list of which mouse buttons are down

	private TransformGroup targetTG;
	private WakeupCondition wakeCond;
	private Topomap3D mainObject;
	private int lastMouseLocX = 0;
	private int lastMouseLocY = 0;
	
	CustomHandleMouse(TransformGroup theTargetTG, Topomap3D theMainObj) {
		this.targetTG = theTargetTG; 
		this.mainObject = theMainObj;
		WakeupCriterion[] wakeCriteria={ new WakeupOnAWTEvent(MouseEvent.MOUSE_PRESSED), 
										 new WakeupOnAWTEvent(MouseEvent.MOUSE_RELEASED),
										 new WakeupOnAWTEvent(MouseEvent.MOUSE_DRAGGED), 
										 new WakeupOnAWTEvent(MouseWheelEvent.MOUSE_WHEEL) };
		wakeCond = new WakeupOr(wakeCriteria);
	}
	
	public void initialize() {
		this.wakeupOn(wakeCond);
	}
	
	// called by Java 3D when appropriate stimulus occurs 
	public void processStimulus(Enumeration criteria) {
		// do what is necessary in response to stimulus
		while (criteria.hasMoreElements()) {
			WakeupOnAWTEvent stimulus = (WakeupOnAWTEvent) criteria.nextElement();
			AWTEvent[] events = stimulus.getAWTEvent();
			for (int i = 0; i < events.length; i++) {
				dealWithMouseEvent((MouseEvent) events[i]);
			}
		}
		this.wakeupOn(wakeCond);
	}
	
	public void dealWithMouseEvent(MouseEvent eve) {
		
		switch (eve.getID()) {
			case MouseEvent.MOUSE_PRESSED:
				lastMouseLocX = eve.getX();
				lastMouseLocY = eve.getY();
				mouseDownMap.set(eve.getButton());
				
				if (mouseDownMap.get(MouseEvent.BUTTON1)) {
					
				} else if (mouseDownMap.get(MouseEvent.BUTTON3)) {
					
				}
				
				break;
			case MouseEvent.MOUSE_RELEASED:
				mouseDownMap.clear(eve.getButton());
				break;
			case MouseEvent.MOUSE_DRAGGED:
				int newMouseLocX = eve.getX();
				int displacementX = newMouseLocX - lastMouseLocX;
				
				int newMouseLocY = eve.getY();
				int displacementY = newMouseLocY - lastMouseLocY;
				
				lastMouseLocX = newMouseLocX;
				lastMouseLocY = newMouseLocY;
				
				mainObject.getCamera().translate(-(float)displacementX / 100.0f, (float)displacementY / 100.0f, 0);
				
				break;
			case MouseWheelEvent.MOUSE_WHEEL:
				// I don't think this is the right way to handle the scroll wheel!
				
				int amountRot = ((MouseWheelEvent)eve).getWheelRotation();
				// Get wheel rotation:
				// Returns the number of "clicks" the mouse wheel was rotated.
				// negative values if the mouse wheel was rotated up (away from) the user,
				// and positive values if the mouse wheel was rotated down (toward) the user
				System.err.println("Mouse wheel was scrolled by " + amountRot + " clicks.");
				// This doesn't seem to work on the Mac, at least not with my mouse.
				break;
			default:
				System.err.println("CustomHandleMouse is handling events that it shouldn't be.");
		}
		
		mainObject.redoView();
		
		//System.err.println("Got a mouse event. " + eve.paramString());
		
	} // end of method dealWithMouseEvent(...)
	
} // end of class CustomHandleMouse 