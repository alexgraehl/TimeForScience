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

// This HUD uses the "NewDawnHUD" Java library (http://www.newdawnsoftware.com)

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.Date;

import javax.media.j3d.*;

import java.awt.*;
import javax.swing.JTextArea;

import org.newdawn.common.hud.*;
import org.newdawn.common.hud.event.*;
import org.newdawn.common.logger.Logger;

public class HUDLegend implements HUDComponentPainter {
    
	private Topomap3D topoObj;
	
    private Graphics2D graphics;
    private HUDComponent component;
    private HUD hud;
	
    private JTextArea textArea;
    
	private Color bgColor   = new Color(0.2f, 0.6f, 0.2f, 0.3f); // Red, green, blue, alpha (0 = transparent, 1 = opaque)
	private Color textColor = new Color(1.0f, 1.0f, 1.0f, 1.0f); // Red, green, blue, alpha (0 = transparent, 1 = opaque)
	
	private static final Color clearColor = new Color(0, 0, 0, 0); // A "clear" (black & transparent) color for erasing with
	
    /** Creates new TestLabel */
    public HUDLegend(HUD hud, Topomap3D theMainObj) {
		
		this.topoObj = theMainObj;
		this.hud = hud;
        
		this.textArea = new JTextArea();
		textArea.setText("Topomap3D\nCurrent settings: (none)");
		
		//bgColor = new Color(Color.BLUE.getRed(), Color.BLUE.getGreen(), Color.BLUE.getBlue(), 255);
		
		float legendRelativeWidth = 1.0f; // 1.0 = full width, 0.0 = none at all
		int   legendHeight = 100; // height in pixels
		MixedSize size = new MixedSize();
		size.addSize(new RelativeSize(legendRelativeWidth, 0f)); // The WIDTH is relative to the size of the window
		size.addSize(new AbsoluteSize(0, legendHeight)); // But the HEIGHT is absolute
        
		MixedPosition position = new MixedPosition();
		position.addPosition(new RelativePosition(0f, 1f)); // put it at the bottom of the screen
		position.addPosition(new AbsolutePosition(0, -legendHeight )); // make sure to move it up so it's visible, though
		
		component = hud.getHUDComponent(1, position, size, this);
		//component.addHUDComponentKeyListener(this); // listen for keypresses (so we can update the HUD)
		//component.requestFocus();
		
		textArea.setBackground(bgColor);
		textArea.setForeground(textColor); // this will be the text color!
		textArea.setLineWrap(true);
		textArea.setWrapStyleWord(true);
    }
    
	// The method that is responsible for painting the the HUDComponent
    public void paint(HUDComponent component) {
		if (component == this.component) {
			this.draw();
		}
    }
	
    public void draw() {
		this.recalculateTextBox();
		
		int width  = component.getSize().getWidth(hud.getWidth());
		int height = component.getSize().getHeight(hud.getHeight());
		
        graphics = (Graphics2D)component.getGraphics();
		
		int fontHeight = graphics.getFontMetrics().getMaxDescent() + graphics.getFontMetrics().getMaxAscent() + graphics.getFontMetrics().getLeading();
		int fontStartX = 2;
		int fontStartY = 2 + graphics.getFontMetrics().getMaxAscent();
		
        graphics.setBackground(clearColor);
        graphics.clearRect(0, 0, width-1, 99);
        //graphics.setColor(bgColor);
        //graphics.fillRect(0, 0, width-1, 99);
        graphics.setColor(Color.BLACK);
        graphics.drawRect(0, 0, width-1, 99);
        graphics.setColor(Color.WHITE);
		textArea.setSize(width-2, height-2);
		textArea.setBorder(new javax.swing.border.BevelBorder(javax.swing.border.BevelBorder.LOWERED));

		int lineCount = textArea.getLineCount();
		int maxLinesInBox = ((int)height - 2) / (int)fontHeight;
		if (lineCount > maxLinesInBox) {
			int startOffset = 0;
			try {
				startOffset = textArea.getLineStartOffset(lineCount - maxLinesInBox - 1);
				String tempText = textArea.getText();
				textArea.setText(tempText.substring(startOffset));
				lineCount = textArea.getLineCount();
				//Logger.log(Logger.DEBUG, "Text area contains " + lineCount + " lines of text");
			} catch (javax.swing.text.BadLocationException exc) {
				System.err.println(exc);
				exc.printStackTrace();
			}
		}
		
		textArea.paint(graphics.create(1,1,width-2, height-2));
        component.imageUpdated();
    }
	
	private void recalculateTextBox() {
		textArea.setText("+/-: Zoom in & out,  T: Hide/show 3D Topomap"  + "\n"
						+"E: Hide/show edges,  N: Hide/show graph nodes"  + "\n"
						+"L: Hide/show graph labels" + "\n"
						+"Mouse: rotate graph (Z: reset camera position)" + "\n"
						+"Y / U / I: Toggle transparency settings" + "\n"
						//+"R: randomize topomap (very very very slow)"
						);
		//+"Vertical Scale: " + topoObj.getTopomapVerticalScale() + "\n"
		//+"Showing graph nodes: "  + topoObj.getShowGraphNodes() + "\n"
		//+"Showing graph edges: "  + topoObj.getShowGraphEdges() + "\n"
		//+"Showing graph labels: " + topoObj.getShowGraphLabels() + "\n"
		//+"Showing topomap: " + topoObj.getShowTopomap() + "\n");
	}
}
