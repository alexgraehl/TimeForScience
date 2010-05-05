;;; show-wspace.el --- Highlight whitespace of various kinds.
;;
;; Filename: show-wspace.el
;; Description: Highlight whitespace of various kinds.
;; Author: Peter Steiner <unistein@isbe.ch>, Drew Adams
;; Maintainer: Drew Adams
;; Copyright (C) 2000-2006, Drew Adams, all rights reserved.
;; Created: Wed Jun 21 08:54:53 2000
;; Version: 21.0
;; Last-Updated: Fri Jul 07 16:48:43 2006 (-25200 Pacific Daylight Time)
;;           By: dradams
;;     Update #: 219
;; URL: http://www.emacswiki.org/cgi-bin/wiki/show-wspace.el
;; Keywords: highlight, whitespace
;; Compatibility: GNU Emacs 20.x, GNU Emacs 21.x, GNU Emacs 22.x
;;
;; Features that might be required by this library:
;;
;;   None
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;; Commentary:
;;
;;    Highlight whitespace of various kinds.
;;
;; To use this library:
;;
;;    Add this to your initialization file (~/.emacs or ~/_emacs):
;;
;;      (require 'show-wspace) ; Load this library.
;;
;; Then you can use commands `toggle-*' (see below) to turn the
;; various kinds of whitespace highlighting on and off in Font-Lock
;; mode.
;;
;; If you want to always use a particular kind of whitespace
;; highlighting, by default, then add the corresponding `highlight-*'
;; command (see below) to the hook `font-lock-mode-hook'.  Then,
;; whenever Font-Lock mode is turned on, so will the whitespace
;; highlighting.
;;
;; For example, you can turn on tab highlighting by default by adding
;; command `highlight-tabs' to `font-lock-mode-hook' in your .emacs
;; file, as follows:
;;
;;     (add-hook 'font-lock-mode-hook 'highlight-tabs)
;;
;;
;; Faces defined here:
;;
;;    `pesche-hardspace', `pesche-space', `pesche-tab'.
;;
;; Commands defined here:
;;
;;    `toggle-hardspace-font-lock', `toggle-tabs-font-lock',
;;    `toggle-trailing-whitespace-font-lock'.
;;
;; Non-interactive functions defined here:
;;
;;    `highlight-hard-spaces', `highlight-tabs',
;;    `highlight-trailing-whitespace'.
;;
;; Internal variables defined here:
;;
;;    `highlight-hard-spaces-p', `highlight-tabs-p',
;;    `highlight-trailing-whitespace-p'.
;;
;; Drew Adams wrote the `toggle-*' commands and `*-p' variables.
;;
;; Peter Steiner wrote the original code that did the equivalent of
;; the `highlight-*' commands here in his `hilite-trail.el'.  The
;; names "pesche" are his.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;; Change log:
;;
;; 2006/04/06 dadams
;;     highlight-*: Use font-lock-add-keywords.  Thanks to Karl Chen.
;; 2006/02/20 dadams
;;     Mentioned in Commentary how to use non-interactively.
;; 2006/01/07 dadams
;;     Added :link for sending bug report.
;; 2006/01/06 dadams
;;     Added defgroup and use it.
;; 2005/12/30 dadams
;;     Removed require of def-face-const.el.
;;     Renamed faces, without "-face".
;; 2005/01/25 dadams
;;     Removed ###autoload for defvars.
;; 2004/06/10 dadams
;;     Fixed minor bug in highlight-* functions.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program; see the file COPYING.  If not, write to
;; the Free Software Foundation, Inc., 51 Franklin Street, Fifth
;; Floor, Boston, MA 02110-1301, USA.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;; Code:

(and (< emacs-major-version 20) (eval-when-compile (require 'cl))) ;; when, push

;;;;;;;;;;;;;;;;;;;;;;;;;

(defgroup Show-Whitespace nil
  "Highlight whitespace of various kinds."
  :group 'convenience :group 'matching
  :link `(url-link :tag "Send Bug Report"
          ,(concat "mailto:" "drew.adams" "@" "oracle" ".com?subject=\
show-wspace.el bug: \
&body=Describe bug here, starting with `emacs -q'.  \
Don't forget to mention your Emacs and library versions."))
  :link '(url-link :tag "Other Libraries by Drew"
          "http://www.emacswiki.org/cgi-bin/wiki/DrewsElispLibraries")
  :link '(url-link :tag "Download"
          "http://www.emacswiki.org/cgi-bin/wiki/show-wspace.el")
  :link '(url-link :tag "Description"
          "http://www.emacswiki.org/cgi-bin/wiki/ShowWhiteSpace#ShowWspace")
  :link '(emacs-commentary-link :tag "Commentary" "show-wspace")
  )

(defface pesche-tab '((t (:background "black" :bold t :underline t :foreground "black")))
  "*Face for highlighting tab characters (`C-i') in Font-Lock mode."
  :group 'Show-Whitespace :group 'font-lock :group 'faces)

(defface pesche-space '((t (:background "red" :bold t :underline t :foreground "black")))
  "*Face for highlighting whitespace at line ends in Font-Lock mode."
  :group 'Show-Whitespace :group 'font-lock :group 'faces)

(defface pesche-hardspace '((t (:background "Green" :underline t)))
  "*Face for highlighting hard spaces (`\040')in Font-Lock mode."
  :group 'Show-Whitespace :group 'font-lock :group 'faces)


(defvar highlight-tabs-p nil
  "Non-nil means font-Lock mode is highlighting TABs (`C-i').")

(defvar highlight-trailing-whitespace-p nil
  "Non-nil means font-Lock mode is highlighting whitespace at line ends.")

(defvar highlight-hard-spaces-p nil
  "Non-nil means font-Lock mode is highlighting hard spaces (`\040').")

;;;###autoload
(defun toggle-tabs-font-lock ()
  "Toggle highlighting of TABs, using face `pesche-tab'."
  (interactive)
  (if highlight-tabs-p
      (remove-hook 'font-lock-mode-hook 'highlight-tabs)
    (add-hook 'font-lock-mode-hook 'highlight-tabs))
  (setq highlight-tabs-p (not highlight-tabs-p))
  (font-lock-mode)(font-lock-mode)
  (message "TAB highlighting is now %s." (if highlight-tabs-p "ON" "OFF")))

;;;###autoload
(defun toggle-hardspace-font-lock ()
  "Toggle highlighting of hard SPACE characters.
Uses face `pesche-hardspace'."
  (interactive)
  (if highlight-hard-spaces-p
      (remove-hook 'font-lock-mode-hook 'highlight-hard-spaces)
    (add-hook 'font-lock-mode-hook 'highlight-hard-spaces))
  (setq highlight-hard-spaces-p (not highlight-hard-spaces-p))
  (font-lock-mode)(font-lock-mode)
  (message "Hard space highlighting is now %s."
           (if highlight-hard-spaces-p "ON" "OFF")))

;;;###autoload
(defun toggle-trailing-whitespace-font-lock ()
  "Toggle highlighting of trailing whitespace.
Uses face `pesche-space'."
  (interactive)
  (if highlight-trailing-whitespace-p
      (remove-hook 'font-lock-mode-hook 'highlight-trailing-whitespace)
    (add-hook 'font-lock-mode-hook 'highlight-trailing-whitespace))
  (setq highlight-trailing-whitespace-p (not highlight-trailing-whitespace-p))
  (font-lock-mode)(font-lock-mode)
  (message "Trailing whitespace highlighting is now %s."
           (if highlight-trailing-whitespace-p "ON" "OFF")))

(defun highlight-tabs ()
  "Highlight tab characters (`C-i')."
  (font-lock-add-keywords nil '(("[\t]+" (0 'pesche-tab t)))))
(defun highlight-hard-spaces ()
  "Highlight hard-space characters (`\040')."
  (font-lock-add-keywords nil '(("[\240]+" (0 'pesche-hardspace t)))))
(defun highlight-trailing-whitespace ()
  "Highlight whitespace characters at line ends."
  (font-lock-add-keywords nil '(("[\040\t]+$" (0 'pesche-space t)))))

;;;;;;;;;;;;;;;;;;;;;;;

(provide 'show-wspace)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; show-wspace.el ends here
