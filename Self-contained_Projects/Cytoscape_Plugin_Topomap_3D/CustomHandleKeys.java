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

import java.util.Enumeration;
import java.util.BitSet;
import java.util.HashMap;

public class CustomHandleKeys extends Behavior {

	private BitSet keyDownMap = new BitSet(); // a list of which keys are up and which are down
	
	private TransformGroup targetTG;
	private WakeupCondition wakeCond;
	private Topomap3D mainObject;
	
	// create CustomHandleKeys - set TG object of change 
	CustomHandleKeys(TransformGroup theTargetTG, Topomap3D theMainObj) { 
		this.targetTG = theTargetTG; 
		this.mainObject = theMainObj;
		
		WakeupCriterion[] wakeCriteria={ new WakeupOnAWTEvent(KeyEvent.KEY_PRESSED), 
										 new WakeupOnAWTEvent(KeyEvent.KEY_RELEASED) };
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
				dealWithKeyEvent((KeyEvent) events[i]);
			}
		}
		this.wakeupOn(wakeCond);
	}
	
	public void handleTimeUpdate() {
		// This deals with keys that have an effect as they are held down.
		// So instead of being per-keypress, they are per-timeunit.
		// Keyboard control of the camera is a good example.
		// Things like toggling of options with the keyboard are NOT handled here,
		// because you want to toggle them once per keypress, not once per 3d drawing frame.
		
		if (keyDownMap.get(61) // +
			|| keyDownMap.get(107)) { // "+" on keypad
			// ZOOM IN
			//mainObject.getCamera().simpleMove(-0.5);
			mainObject.getCamera().translate(0, 0, -0.2f);
		}
		if (keyDownMap.get(45) // -
			|| keyDownMap.get(109)) { // "-" on keypad
			// ZOOM OUT
			mainObject.getCamera().translate(0, 0, +0.2f);
		}
		
		if (keyDownMap.get(37) // left arrow
			|| keyDownMap.get(100)) { // "4" on numeric keypad
			Transform3D temp = new Transform3D();
			temp.rotZ(-0.1f);
			Transform3D theOldTransform = new Transform3D();
			this.targetTG.getTransform(theOldTransform); // copies "targetTG"'s transform into theGroupTransform
			temp.mul(theOldTransform); // modify temp
			this.targetTG.setTransform(temp); // figures out the new transformation matrix
		}
		if (keyDownMap.get(38) // up arrow
			|| keyDownMap.get(104)) { // "8" on numeric keypad
			Transform3D temp = new Transform3D();
			temp.rotY(0.1f);
			Transform3D theOldTransform = new Transform3D();
			this.targetTG.getTransform(theOldTransform); // copies "targetTG"'s transform into theGroupTransform
			temp.mul(theOldTransform); // modify temp
			this.targetTG.setTransform(temp); // figures out the new transformation matrix
		}
		if (keyDownMap.get(39) // right arrow
			|| keyDownMap.get(102)) { // "6" on numeric keypad
			Transform3D temp = new Transform3D();
			temp.rotZ(0.1f);
			Transform3D theOldTransform = new Transform3D();
			this.targetTG.getTransform(theOldTransform); // copies "targetTG"'s transform into theGroupTransform
			temp.mul(theOldTransform); // modify temp
			this.targetTG.setTransform(temp); // figures out the new transformation matrix
		}
		if (keyDownMap.get(40) // down arrow
			|| keyDownMap.get(98)) { // "2" on numeric keypad
			Transform3D temp = new Transform3D();
			temp.rotY(-0.1f);
			Transform3D theOldTransform = new Transform3D();
			this.targetTG.getTransform(theOldTransform); // copies "targetTG"'s transform into theGroupTransform
			temp.mul(theOldTransform); // modify temp
			this.targetTG.setTransform(temp); // figures out the new transformation matrix
		}
	}
	
	public void dealWithKeyEvent(KeyEvent theEvent) {
		int eventID = theEvent.getID();
		if (eventID == KeyEvent.KEY_RELEASED) {
			int keyCode = ((KeyEvent)theEvent).getKeyCode();
			keyDownMap.clear(keyCode); // update the keyDownMap (key is no longer down)
			
		} else if (eventID == KeyEvent.KEY_PRESSED) {
			char keyChar = ((KeyEvent)theEvent).getKeyChar();
			int keyLocationOnKbd  = ((KeyEvent)theEvent).getKeyLocation();
			int keyCode = ((KeyEvent)theEvent).getKeyCode();
			
			keyDownMap.set(keyCode); // update the keyDownMap
			
			switch (keyChar) {
				case '?':  case '/':
					System.err.println("Should show help");
					break;
					
				case 'r':  case 'R':
					//javax.media.j3d.Locale myLoc = mainObject.getUniverse().getLocale();
					//mainObject.testMark(mainObject, 1, 1, 0, 3);
					//mainObject.redoView();
					break;
					
				case 'z':  case 'Z':
					mainObject.getCamera().moveCameraTo(new Point3d(0, 0, 2),
														new Point3d(0, 0, 0));
					mainObject.redoView();
					break;
				
				case 'n':  case 'N': // Toggle showing of graph nodes
					mainObject.setShowGraphNodes(!mainObject.getShowGraphNodes());
					mainObject.redoBranches();
					//System.err.println("Toggled 'n' to " + mainObject.getShowGraphNodes());
					break;
				case 'e':  case 'E': // Toggle showing of graph edges
					mainObject.setShowGraphEdges(!mainObject.getShowGraphEdges());
					mainObject.redoBranches();
					//System.err.println("Toggled 'e' to " + mainObject.getShowGraphEdges());
					break;
				case 'l':  case 'L': // Toggle showing of graph labels
					mainObject.setShowGraphLabels(!mainObject.getShowGraphLabels());
					mainObject.redoBranches();
					//System.err.println("Toggled 'l' to " + mainObject.getShowGraphLabels());
					break;
				case 't':  case 'T': // Toggle showing of topomap
					mainObject.setShowTopomap(!mainObject.getShowTopomap());
					mainObject.redoBranches();
					//System.err.println("Toggled 't' to " + mainObject.getShowTopomap());
					break;
				case 'y': case 'Y': // Toggle transparency settings (screen door)
					mainObject.getTerrainShape3D().setAppearance(LandscapeShape3D.makeLandscapeAppearance(javax.media.j3d.TransparencyAttributes.SCREEN_DOOR, 0.5f));
					mainObject.redoView();
					break;
				case 'u': case 'U': // Toggle transparency settings
					mainObject.getTerrainShape3D().setAppearance(LandscapeShape3D.makeLandscapeAppearance(javax.media.j3d.TransparencyAttributes.NICEST, 0.5f));
					mainObject.redoView();
					break;
				case 'i': case 'I': // Toggle transparency settings (opaque)
					mainObject.getTerrainShape3D().setAppearance(LandscapeShape3D.makeLandscapeAppearance(javax.media.j3d.TransparencyAttributes.NONE, 0.0f));
					mainObject.redoView();
					break;
			}
			
			System.out.println("Key pressed: char: " + keyChar + " loc: " + keyLocationOnKbd + " keycode: " + keyCode);
			
		} else {
			System.out.println("CustomHandleKeys is handling events that it shouldn't be.");
		}
		
	} // end of method dealWithKeyEvent(...)
	
} // end of class CustomHandleKeys 