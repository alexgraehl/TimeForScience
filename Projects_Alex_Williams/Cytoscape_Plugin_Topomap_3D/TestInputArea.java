//package org.newdawn.common.hud.demo;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.Date;

import javax.swing.JTextArea;

import org.newdawn.common.hud.*;

import org.newdawn.common.hud.event.FocusListener;
import org.newdawn.common.hud.event.KeyListener;
import org.newdawn.common.hud.event.MouseListener;
import org.newdawn.common.hud.AbsolutePosition;
import org.newdawn.common.hud.AbsoluteSize;
import org.newdawn.common.hud.HUD;
import org.newdawn.common.hud.HUDComponent;
import org.newdawn.common.hud.HUDComponentListener;
import org.newdawn.common.hud.HUDComponentPainter;
import org.newdawn.common.hud.MixedPosition;
import org.newdawn.common.hud.MixedSize;
import org.newdawn.common.hud.RelativePosition;
import org.newdawn.common.hud.RelativeSize;
import org.newdawn.common.logger.Logger;

/**
 *
 * @author  Jeremy
 * @version 
 */
public class TestInputArea implements MouseListener, KeyListener, FocusListener, HUDComponentPainter {
    
    private Color backGroundColour;
    private Color backGroundColourFaded;
    private Color backGroundColourAlphad;
    private Color clearColour;
    private Graphics2D graphics;
    private HUDComponent component;
    private HUD hud;
    private TestDisplayArea displayArea;
    
    private JTextArea textArea = new JTextArea();
    
    /** Creates new TestLabel */
    public TestInputArea(HUD hud, TestDisplayArea displayArea) {
	
	this.hud = hud;
	this.displayArea = displayArea;
        
        backGroundColour = new Color(Color.BLUE.getRed(), Color.BLUE.getGreen(), Color.BLUE.getBlue(), 255);
        backGroundColourFaded = new Color(Color.BLUE.getRed(), Color.BLUE.getGreen(), Color.BLUE.getBlue(), 64);
        backGroundColourAlphad = new Color(Color.BLUE.getRed(), Color.BLUE.getGreen(), Color.BLUE.getBlue(), 0);
/*        backGroundColour = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 255);
        backGroundColourFaded = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 128);
        backGroundColourAlphad = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 0);*/
        clearColour = new Color(Color.BLACK.getRed(), Color.BLACK.getGreen(), Color.BLACK.getBlue(), 0);
        MixedPosition position = new MixedPosition();
        
        position.addPosition(new RelativePosition(0f,1f));
        position.addPosition(new AbsolutePosition(0,-30));
	
	MixedSize size = new MixedSize();
	
	size.addSize(new RelativeSize(0.5f,0f));
	size.addSize(new AbsoluteSize(0,31));
        
        component = hud.getHUDComponent(1, position, size, this);
        
        //component.setTransparency(0f);
        
	component.addHUDComponentMouseListener(this);
	component.addHUDComponentKeyListener(this);
	component.addHUDComponentFocusListener(this);
	
	component.requestFocus();
	
	textArea.setBackground(backGroundColourFaded);
	textArea.setForeground(Color.WHITE);
	textArea.setLineWrap(true);
	textArea.setWrapStyleWord(true);
        
    }
    
    private void draw() {
        
	int width = component.getSize().getWidth(hud.getWidth());
	int height = component.getSize().getHeight(hud.getHeight());
	
	//Logger.log(Logger.DEBUG, "TestInputArea width: " + width + " height: " + height);
	
        graphics = (Graphics2D)component.getGraphics();
	
	int fontHeight = graphics.getFontMetrics().getMaxDescent() + graphics.getFontMetrics().getMaxAscent() + graphics.getFontMetrics().getLeading();
	int fontStartX = 2;
	int fontStartY = 2 + graphics.getFontMetrics().getMaxAscent();

        graphics.setBackground(clearColour);
        
        graphics.clearRect(0, 0, width-1, 99);
        
        //graphics.setColor(backGroundColourFaded);

        //graphics.fillRect(0, 0, width-1, 99);

        graphics.setColor(Color.BLACK);
        
        graphics.drawRect(0, 0, width-1, 99);
        
        graphics.setColor(Color.WHITE);
	
	textArea.setSize(width-2, height-2);
	
	textArea.setBorder(new javax.swing.border.BevelBorder(javax.swing.border.BevelBorder.LOWERED));
	
	int lineCount = textArea.getLineCount();
	if(lineCount>1) {
	    int startOffset = 0;
	    try {
		startOffset = textArea.getLineStartOffset(lineCount - 1);
		String tempText = textArea.getText();
		textArea.setText(tempText.substring(startOffset));
		lineCount = textArea.getLineCount();
		//Logger.log(Logger.DEBUG, "Text area contains " + lineCount + " lines of text");
	    } catch (javax.swing.text.BadLocationException e) {
		Logger.log(Logger.WARNING, "Hmmm", e);
	    }
	}
	
	textArea.paint(graphics.create(1,1,width-2, height-2));
        
        component.imageUpdated();
    }

    /** Notification that the mouse has been clicked inside this component
     *
     * @param component The component selected
     * @param x The x coord where the mouse was clicked
     * @param y The y coord where the mouse was clicked
     */
    public void mouseClicked(org.newdawn.common.hud.event.MouseEvent mouseEvent) {
	component.requestFocus();
	mouseEvent.consume();
    }
    
    /**
     * Notification that a key has been pressed
     *
     * @param keyEvent The key event for this key type
     */
    public void keyTyped(org.newdawn.common.hud.event.KeyEvent keyEvent) {
	String origString = textArea.getText();
	char keyChar = keyEvent.getKeyChar();
	String newString = origString + keyChar;
	
	if(keyChar == '\b') {
	    if(origString.length()>0) {
		int length = origString.length();
		newString = origString.substring(0, length - 1);
		textArea.setText(newString);
		draw();
		keyEvent.consume();
		return;
	    }
	}
	if(keyChar == '\n') {	    
	    if(origString.length()>0) {
		displayArea.appendText(origString + "\n");
		textArea.setText("");
		draw();
	    }
	    component.surrenderFocus();
	    keyEvent.consume();
	    return;
	}
	
	textArea.setText(newString);
	draw();
	keyEvent.consume();
    }
    
    /** Notification that focus has been gained
     */
    public void focusGained() {
	hud.getEventManager().setHUDExclusiveMode(true);
    }
    
    /** Notification that focus has been lost
     */
    public void focusLost() {
	hud.getEventManager().setHUDExclusiveMode(false);
    }
    
    /**
     * Notification that a key has been pressed
     *
     * @param keyEvent The key event for this key type
     */
    public void keyPressed(org.newdawn.common.hud.event.KeyEvent keyEvent) {
    }
    
    /**
     * Notification that a key has been pressed
     *
     * @param keyEvent The key event for this key type
     */
    public void keyReleased(org.newdawn.common.hud.event.KeyEvent keyEvent) {
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
