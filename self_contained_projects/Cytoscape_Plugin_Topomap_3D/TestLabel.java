//package org.newdawn.common.hud.demo;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.Date;

import org.newdawn.common.hud.AbsolutePosition;
import org.newdawn.common.hud.AbsoluteSize;
import org.newdawn.common.hud.event.MouseListener;
import org.newdawn.common.hud.HUD;
import org.newdawn.common.hud.HUDComponent;
import org.newdawn.common.hud.HUDComponentListener;
import org.newdawn.common.hud.HUDComponentPainter;
import org.newdawn.common.hud.MixedPosition;
import org.newdawn.common.hud.RelativePosition;
import org.newdawn.common.logger.Logger;

/**
 *
 * @author  Jeremy
 * @version 
 */
public class TestLabel implements HUDComponentPainter, MouseListener {
    
    private Color backGroundColour;
    private Color backGroundColourFaded;
    private Color backGroundColourAlphad;
    private Color clearColour;
    private Graphics2D graphics;
    private HUDComponent component;

    /** Creates new TestLabel */
    public TestLabel(HUD hud) {
        
        backGroundColour = new Color(Color.BLUE.getRed(), Color.BLUE.getGreen(), Color.BLUE.getBlue(), 255);
        backGroundColourFaded = new Color(Color.BLUE.getRed(), Color.BLUE.getGreen(), Color.BLUE.getBlue(), 128);
        backGroundColourAlphad = new Color(Color.BLUE.getRed(), Color.BLUE.getGreen(), Color.BLUE.getBlue(), 0);
/*        backGroundColour = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 255);
        backGroundColourFaded = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 128);
        backGroundColourAlphad = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 0);*/
        clearColour = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 0);
        MixedPosition position = new MixedPosition();
        
        position.addPosition(new RelativePosition(1f,1f));
        position.addPosition(new AbsolutePosition(-250,-100));
        
        component = hud.getHUDComponent(1, position, new AbsoluteSize(250,100), this);
        
        //component.setTransparency(0f);
        
	component.addHUDComponentMouseListener(this);
        
    }
    
    private void draw() {
        
        graphics = (Graphics2D)component.getGraphics();

        graphics.setBackground(clearColour);
        
        graphics.clearRect(0,0,249,99);
        
        graphics.setColor(backGroundColour);

        graphics.fillRect(0,0,90,99);

        graphics.setColor(backGroundColourFaded);

        graphics.fillRect(90,0,90,99);

        graphics.setColor(backGroundColourAlphad);

        graphics.fillRect(180,0,69,99);

        graphics.setColor(Color.WHITE);
        
        graphics.drawRect(0,0,249,99);
        
        graphics.setColor(Color.WHITE);
        
        graphics.drawString("Test: " + new Date(), 2 ,2 + graphics.getFontMetrics().getAscent() );
        
        component.imageUpdated();
    }

    /** Notification that the mouse has been clicked inside this component
     *
     * @param component The component selected
     * @param x The x coord where the mouse was clicked
     * @param y The y coord where the mouse was clicked
     */
    public void mouseClicked(org.newdawn.common.hud.event.MouseEvent mouseEvent) {
        Logger.log(Logger.DEBUG, "Mouse clicked at " + mouseEvent.getXLocation() + ":" + mouseEvent.getYLocation() + " and consumed in TestLabel");
	component.requestFocus();
	mouseEvent.consume();
    }
    
    /** The method that is responsible for painting the the HUDComponent
     *
     * @param component The component that needs painting too
     */
    public void paint(HUDComponent component) {
	if(component == this.component) {
	    draw();
	}
    }
    
}
