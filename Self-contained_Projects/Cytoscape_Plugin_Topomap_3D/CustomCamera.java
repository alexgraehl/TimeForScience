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


import javax.media.j3d.*;
import javax.vecmath.*;

public class CustomCamera {
	// The three "orient" points/vectors should all be unit length (i.e., normalized).
	public Point3d  orientEyePt;	// Note that the "eye point" is NOT actually where the camera is located!
	public Point3d  orientLookAtPt;
	public Vector3d orientUpVec; // <-- which way is "up" from the camera's perspective
	
	public Point3d location; // <-- this is the ACTUAL location of the camera
	
	public Transform3D translationT3D; // location of camera
	public Transform3D orientationT3D; // ... and where it's looking
	
	private Transform3D tempT3D = new Transform3D();
	
	CustomCamera() {
		orientEyePt    = new Point3d();
		orientLookAtPt = new Point3d();
		orientUpVec    = new Vector3d();
		location       = new Point3d();
		
		translationT3D = new Transform3D();
		orientationT3D = new Transform3D();
		
	}
	
	//private float getPitch() {
		// Pitch is the up-down tilt of the "nose" of the camera.
		//Point2f groundDistPt = new Point2f(orientLookAtPt.x, orientLookAtPt.y);
		//groundDistPt.sub(new Point2f(orientEyePt.x, orientEyePt.y);
		//float groundDistance = groundDistPt.normalize();
		
		//return
	//}
	
	public Transform3D getTotalTransformation() {
		tempT3D.mul(this.translationT3D, this.orientationT3D);
		return tempT3D;
	}
	
	public void setOrientation(Point3d eyePt, Point3d lookAtPt, Vector3d upVec) {
		orientEyePt    = eyePt;
		orientLookAtPt = lookAtPt;
		orientUpVec    = upVec;
		this.orientationT3D.lookAt(eyePt, lookAtPt, upVec);
	}
	
	public void setLocation(Point3d theLocPt) {
		this.location = theLocPt;
		translationT3D.set(new Vector3d(theLocPt));
	}
	
	public void simpleMove(double newDist) {
		// Moves the point in/out along the angle of view
		
		Point3d newLocation = new Point3d();
		newLocation.add(getNormalizedViewDirection(orientEyePt, orientLookAtPt),
						this.location);
		setLocation(this.location);
	}
	
	public void translate(float x, float y, float z) {
		this.setLocation(new Point3d(this.location.x + x, 
									 this.location.y + y, 
									 this.location.z + z));
	}
	
	public void moveCameraTo(Point3d cameraPoint, Point3d lookAtPoint) {
		// The camera is moved to the cameraPoint, and it is set so that it's looking at the lookAtPoint
		// The up vector is set in a way that I haven't figured out yet
		
		this.setLocation(cameraPoint); // easy to move it to the right place...
		
		// The lookingDirectionVec must be NORMALIZED.
		// You can't just put the eyePt and look-at point in normally.
		// Instead you have to figure out which way the eye would be looking in order
		// to look at that point, and put that in instead.
		// Note that the "eye point" is NOT actually where the camera is located!
		this.setOrientation(new Point3d(0, 0, 0), // eye point
							new Point3d( getNormalizedViewDirection(cameraPoint, lookAtPoint) ), // looking at point
							new Vector3d(0, 1, 0)); // up vector
		// The up vector here is not actually smart at all. sometimes it even makes no sense
		// but, this can be fixed later.
	}
	
	public Vector3d getNormalizedViewDirection(Point3d from, Point3d to) {
		Vector3d lookingDirectionVec = new Vector3d();
		lookingDirectionVec.sub(to, from);
		lookingDirectionVec.normalize();
		return lookingDirectionVec;
	}
}
