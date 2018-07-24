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

import java.awt.*;
import java.awt.event.*;
import javax.media.j3d.*;
import java.util.Enumeration;

public class CustomTimer extends Behavior {
	// This class is just a timer that comes up every so often (as specified by numMilliseconds in the constructor) and calls an update method in the original object. Could be used to animate things, or move the user X amount for every second a key is held down (instead of X amount for every press)
	
	private int internalCounter = 0;
	
	private Topomap3D mainObject;
	private WakeupCondition wakeCond;
	
	CustomTimer(long numMilliseconds, Topomap3D theMainObj) { 
		this.mainObject = theMainObj;
		
		WakeupCriterion[] wakeCriteria={ new WakeupOnElapsedTime((long)numMilliseconds),
										 new WakeupOnElapsedFrames((int)numMilliseconds)};	
		
		// we use WakeupOnElapsedFrames as well, because WakeupOnElapsedTime doesn't
		// seem to always work. It's just a backup to ensure that the callback actually happens.
		// I don't know why, but just the WakeupOnElapsedTime condition fails after about 10-60 seconds
		// of playing around with the keyboard controls (processStimulus is no longer called after that point)
		wakeCond = new WakeupOr(wakeCriteria);
	}
	
	public void initialize() {
		this.wakeupOn(wakeCond);
	}
	
	// called by Java 3D when appropriate stimulus occurs 
	public void processStimulus(Enumeration criteria) {
		internalCounter++;
		//if (internalCounter % 10 == 0) {
		//	System.out.println("Internal counter: " + internalCounter);
		//}
		this.mainObject.handleTimeUpdate();
		this.wakeupOn(wakeCond);
	}
	
}