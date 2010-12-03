//package org.newdawn.common.hud.demo;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.Date;

import org.newdawn.common.hud.AbsolutePosition;
import org.newdawn.common.hud.event.MouseListener;
import org.newdawn.common.hud.HUD;
import org.newdawn.common.hud.HUDComponent;
import org.newdawn.common.hud.HUDComponentListener;
import org.newdawn.common.hud.HUDComponentPainter;
import org.newdawn.common.hud.RelativeSize;
import org.newdawn.common.logger.Logger;

/**
 *
 * @author  Jeremy
 */
public class TestGrid implements MouseListener, HUDComponentPainter {
    
    private Color backGroundColourFaded = new Color(Color.GRAY.getRed(), Color.GRAY.getGreen(), Color.GRAY.getBlue(), 128);
    private Color clearColour = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 0);
    HUDComponent component;
    HUD hud;

    /** Creates new TestLabel */
    public TestGrid(HUD hud) {
        this.hud = hud;
        
        component = hud.getHUDComponent(0,new AbsolutePosition(0,0) , new RelativeSize(1f,1f), this);
	component.addHUDComponentMouseListener(this);
        
        drawComponent();
    }
        
    private void drawComponent() {
        
        int width = component.getSize().getWidth(hud.getWidth());
        int height = component.getSize().getHeight(hud.getHeight());

        //System.out.println(" --->>> Drawing to component<<<---");
        
        //component.setTransparency(0f);
        
        Graphics2D graphics = (Graphics2D)component.getGraphics();
        
        graphics.setBackground(clearColour);
        
        graphics.clearRect(0,0,width,height);
        
        graphics.setColor(clearColour);

        graphics.fillRect(0,0,width,height);

        graphics.setColor(Color.RED);
        
        //graphics.drawRect(0,0, width -1, height -1);
        
        //System.out.println("--->>> Drawing grids <<<---");
        for (int i=0; ((10*i) < (width/2)) ; i++) {
            int x = (i*10);
            int y = (i*10);
            int drawWidth = ((width) -1 - ((i*2)*10));
            int drawHeight = ((height) -1 - ((i*2)*10));
            //System.out.println("--->>> Drawing grid <<<---" + x + " " + y + " " + drawWidth + " " + drawHeight);
            graphics.drawRect(x, y, drawWidth, drawHeight);
        }

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
        Logger.log(Logger.DEBUG, "Mouse clicked at " + mouseEvent.getXLocation() + ":" + mouseEvent.getYLocation() + " in TestGrid");
	Logger.log(Logger.DEBUG, "Mouse event: " + mouseEvent.toString());
    }
    
    /** The method that is responsible for painting the the HUDComponent
     *
     * @param component The component that needs painting too
     */
    public void paint(HUDComponent component) {
	if(component == this.component) {
	    drawComponent();
	}
    }
    
}
