// Topomap3D: A 3-D graph-viewer Cytoscape plug-in
// This source file is part of the Topomap3D plug-in.

// Alex Williams
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

import javax.media.j3d.Group;
import javax.media.j3d.BranchGroup;

public class ReattachableBranchGroup extends javax.media.j3d.BranchGroup {
	// This class is almost the same as a regular Java 3D BranchGroup, but it remembers
	// its last parent, so you can easily reattach it.

	private Group lastParent = null; // pointer to the last group that was a parent of this BranchGroup
	private Group currentParent = null; // required, because you can't ask this group "getParent()" when the branch is LIVE (so you have to store it yourself)
	
	ReattachableBranchGroup() {
		super();
		this.setCapability(BranchGroup.ALLOW_DETACH);
		this.setCapability(Group.ALLOW_CHILDREN_READ);
		this.setCapability(Group.ALLOW_CHILDREN_WRITE);
		this.setCapability(Group.ALLOW_CHILDREN_EXTEND);
	}
	
	public void addToParent(Group addOurselvesToThisParent) {
		addOurselvesToThisParent.addChild(this);
		this.currentParent = addOurselvesToThisParent;
		lastParent = addOurselvesToThisParent;
	}
	
	public boolean isAttached() {
		//System.err.println("Attachment status: " + (currentParent == null));
		return (this.currentParent != null); // Does this node have a parent?
	}
	
	public void detach() {
		if (this.isAttached()) {
			lastParent = currentParent; // save the last parent
			currentParent = null; // no current parent
			super.detach(); // call the parent class's detach method in order to ACTUALLY detach this branch
		} else {
			// no parent, so don't deatch from anything
		}
	}
	
	public void reattach() {
		if (this.isAttached()) {
			// Already attached! So don't do anything...
			//System.err.println("Trying to re-attach a ReattachableBranchGroup that is already attached.");			
		} else {
			this.addToParent(lastParent); // add this group back to the last parent...
		}
	}
	
	public void setAttach(boolean shouldBeAttached) {
		// Lets you call whether to attach or detach this graph based on a boolean argument.
		if (shouldBeAttached) {
			this.reattach();
		} else {
			this.detach();
		}
	}
	
}
