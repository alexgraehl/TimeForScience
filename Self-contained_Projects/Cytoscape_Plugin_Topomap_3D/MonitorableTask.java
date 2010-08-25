
// redistributed Cytoscape code

/* * Redistributed Date: Mar.25.2005
 * * by : Steven Maere
 * * Copyright (c) 2005 Flanders Interuniversitary Institute for Biotechnology (VIB)
 * *
 * * This program is free software; you can redistribute it and/or modify
 * * it under the terms of the GNU General Public License as published by
 * * the Free Software Foundation; either version 2 of the License, or
 * * (at your option) any later version.
 * *
 * * This program is distributed in the hope that it will be useful,
 * * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * * The software and documentation provided hereunder is on an "as is" basis,
 * * and the Flanders Interuniversitary Institute for Biotechnology
 * * has no obligations to provide maintenance, support,
 * * updates, enhancements or modifications.  In no event shall the
 * * Flanders Interuniversitary Institute for Biotechnology
 * * be liable to any party for direct, indirect, special,
 * * incidental or consequential damages, including lost profits, arising
 * * out of the use of this software and its documentation, even if
 * * the Flanders Interuniversitary Institute for Biotechnology
 * * has been advised of the possibility of such damage. See the
 * * GNU General Public License for more details.
 * *
 * * You should have received a copy of the GNU General Public License
 * * along with this program; if not, write to the Free Software
 * * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * */

import java.lang.String;

/**
 * Classes that perform long tasks (like graph algorithms) can implement this interface
 * so that they can be monitored by a GUI like a <code>javax.swing.plaf.ProgressBarUI</code> or a
 * <code>javax.swing.ProgressMonitor</code>
 * Note: this was copied from giny.util because it is being phased out.  Eventually
 * the layout API will be available to use (TODO: remove when layout API is available)
 */
// Ohmagawd!  This is a wonderfully complex and poorly documented interface.
public interface MonitorableTask {

  /**
   * @return <code>true</code> if the task is done, false otherwise
   */
  public boolean isDone();

  /**
   * @return the current progress
   */
  public int getCurrentProgress();

  /**
   * @return the total length of the task
   */
  public int getLengthOfTask();

  /**
   * @return a <code>String</code> describing the task being performed
   */
  public String getTaskDescription();

  /**
   * @return a <code>String</code> status message describing what the task
   * is currently doing (example: "Completed 23% of total.", "Initializing...", etc).
   */
  public String getCurrentStatusMessage ();

  /**
   * Starts doing the task in a separate thread so that the GUI stays responsive
   *
   * @param return_when_done if <code>true</code>, then this method will return only when
   * the task is done, else, it will return immediately after spawning the thread that
   * performs the task
   */
  public void start (boolean return_when_done);

  /**
   * Stops the task if it is currently running.
   */
  public void stop();

  /**
   * @return <code>true</code> if the task was canceled before it was done
   * (for example, by calling <code>MonitorableSwingWorker.stop()</code>,
   * <code>false</code> otherwise
   */
  // TODO: Not sure if needed
  public boolean wasCanceled ();

}//class MonitorableTask
