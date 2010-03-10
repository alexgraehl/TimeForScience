
;; -*-Lisp-*-  <-- Tells emacs what syntax highlighting to use ; Important if the line ending does not indicate the file type.

;; How to set a keybinding interactively:
;; 1. M-x: global-set-key
;; 2. Type the key combination you want
;; 3. Use C-x ESC ESC (repeat-complex-command) to see the apropiate command. Example:
;;     (global-set-key (quote [C-tab]) (quote my-func))

;; M-x apropos: searches command documentation (very limited and unreliable)
;; M-x command-apropos: searches for the literal text of a command

;; http://xahlee.org/emacs/keyboard_shortcuts.html

;;(setq debug-on-error t)

;;(setq tags-file-name "~/cvs-local/TAGS")
;;(visit-tags-table "~/cvs-local/TAGS")

;; http://chumsley.org/download/markdown.el NA ND

;; Emacs regexp: http://www.cs.utah.edu/dept/old/texinfo/emacs18/emacs_17.html

;; List of all emacs colors: http://www.geocities.com/kensanata/colors.html

(defvar ACTIVE_MODELINE       "red")   ;; Background for the active modeline
(defvar ACTIVE_MODELINE_TEXT  "yellow") ;; Text in the currently-active modeline

(defvar INACTIVE_MODELINE       "blue") ;; Background for the inactive modeline
(defvar INACTIVE_MODELINE_TEXT  "white") ;; Text in the inactive modelines

(defvar BUFFER_FILENAME_TEXT "blue")   ;; colors for the filename
(defvar BUFFER_FILENAME_BG   "yellow")

(defvar CURSOR_TEXT  "yellow")
(defvar CURSOR_BG    "red")

(defvar BACKGROUND   "black")

(defvar SEARCH_RESULT_HIGHLIGHT_TEXT "red")
(defvar SEARCH_RESULT_HIGHLIGHT_BG   "yellow")

(defvar SEARCH_RESULT_OTHER_TEXT SEARCH_RESULT_HIGHLIGHT_BG)
(defvar SEARCH_RESULT_OTHER_BG   SEARCH_RESULT_HIGHLIGHT_TEXT)

(defvar COMMENT_FOREGROUND "blue") ;"dark slate gray") ; was dark red before
(defvar COMMENT_BACKGROUND nil)

(defvar STRING_TEXT "green") ;;"black")
(defvar STRING_BG   "blue")  ;;"green")

(defvar HASH_TEXT     "yellow") ; only usable in cperl mode
(defvar ARRAY_TEXT    "yellow") ; only usable in ceperl mode
(defvar VARIABLE_TEXT "cyan")

(defvar VARIABLE_BG   "blue")

(defvar FUNCTION_TEXT "cyan")
(defvar FUNCTION_BG   "blue")

(defvar TYPE_TEXT     "yellow")
(defvar TYPE_BG       nil)       ;"orange4")

(defvar SELECTED_TEXT "yellow") ; Highlighted region / selection
(defvar SELECTED_BG   "red")

(defvar ERR_TEXT "black")
(defvar ERR_BG   "yellow")

(defvar BUILTIN_TEXT "black")
(defvar BUILTIN_BG   "magenta")

(defvar DOC_TEXT "red")
(defvar DOC_BG   STRING_BG)

(defvar MINIBUFFER_TEXT "yellow")
(defvar MINIBUFFER_BG   "red")

(defvar KEYWORD_TEXT "violet")
(defvar KEYWORD_BG   nil)

(defface agwCustomDoubleLineFace
  '((t (:foreground "magenta" :background "black" :underline t))) "To highlight linebreak-style rows of == characters.")

(defface agwLineNumberFace
  '((t (:foreground "black" :background "blue" :inverse-video nil)))
  "Face for the line numbers at the beginning of each line")

(defface agwPositiveNumberFace
  '((t (:foreground "green"))) "To highlight positive numbers")

(defface agwZeroNumberFace
  '((t (:foreground "cyan")))
  "To highlight zeros")

(defface agwNegativeNumberFace
  '((t (:foreground "red")))
  "To highlight negative numbers")


(defface agwPythonPassFace
  '((t (:bold f  :foreground "black"  :background "magenta"  :inverse-video nil)))
  "Face for PASS in Python."  :group 'agwFaces)

(defface agwTempFilenameFace
  '((t (:bold t  :foreground "black"  :background "green"  :inverse-video nil)))
  "Face for temp files in MAKE"  :group 'agwFaces)

(defface agwArrayFace
  '((t ( ;;:inherit font-lock-warning-face
	:bold nil :background nil :foreground "green"))) "Indicates arrays/vectors." :group 'agwFaces)

(defface agwListFace
  '((t ( ;; inherit some-other-face
	:bold nil :foreground "green"))) "Indicates which variables are an R list, based on the name." :group 'agwFaces)

(defface agwAssertionFace
  '((t ( ;; inherit some-other-face
	:bold nil :foreground "yellow" :background "black" :underline nil))
    )
  "Highlights assertions."
  :group 'agwFaces)

(defface agwMakeGlobalVarFace
  '((t ( ;; inherit some-other-face
	:bold nil :foreground "black" :background "magenta" :underline nil))
    )
  "Indicates which variables are GLOBAL variables, assuming they start with 'gv'"
  :group 'agwFaces)



(defface agwIndent1Face
  '((t ( ;; inherit some-other-face
	:bold nil :background "red" :foreground nil :underline nil))
    )  ""  :group 'agwFaces)

(defface agwIndent2Face
  '((t ( ;; inherit some-other-face
	:bold nil :background "blue" :foreground nil :underline nil))
    )  ""  :group 'agwFaces)

(defface agwIndent3Face
  '((t ( ;; inherit some-other-face
	:bold nil :background "red" :foreground nil :underline nil))
    )  ""  :group 'agwFaces)

(defface agwIndent4Face
  '((t ( ;; inherit some-other-face
	:bold nil :background "blue" :foreground nil :underline nil))
    )  ""   :group 'agwFaces)






;;(require 'cl nil t) ;; a rare necessary use of REQUIRE <-- (no idea what this does)
;;(defvar *emacs-load-start* (current-time)) ; <-- uncomment this to figure out how long it took to run this .emacs

;; Get current system type
					;(defun insert-system-type()  (interactive)  "Get current system type" (insert (format "%s" system-type)))
(defun insert-system-type()  (format "%s" system-type))
(defun system-type-is-darwin () (string-equal system-type "darwin"))
(defun system-type-is-gnu () (string-equal system-type "gnu/linux"))

(defun show-system-type() (interactive) (message (insert-system-type)))

(setq should-load-ess nil) ;; <-- initialize a new variable if we haven't loaded ess, then this *remains* "nil"


(setq load-path (cons "~/.emacs.d" load-path))

(if (system-type-is-darwin)
    (progn
      (setq load-path (cons "/usr/local/share/emacs/site-lisp" load-path))
      (setq should-load-ess t)
      (require 'ess-site nil t) ;; <-- this is super slow!!!
      ))

(if (system-type-is-gnu)
    (progn
      (setq load-path (cons "/usr/local/share/emacs/site-lisp" load-path))
      (setq should-load-ess t)
      (require 'ess-site nil t) ;; <-- this can be super slow!!!
      (require 'show-wspace nil t)			; Show whitespace!
      ))

;;(autoload 'ruby-mode "ruby-mode" "Major mode for editing ruby scripts." t)
;;(load (expand-file-name "~/.emacs.d/cperl-mode.el.6.2")) ; <-- note: compiled version!
;;(autoload 'cperl-mode "cperl-mode" "Fancier mode Perl highlighting." t)

;; ##########################################################################
;; ESS (Emacs Speaks Statistics == "R-project" R.app r-project.org R stuff)
;;(if (system-type-is-gnu) ; load only on the GNU filesystem at work
;;    (load (expand-file-name "~/.emacs.d/ess/lisp/ess-site")))
;; Parameters to load that will make it not complain/warn: 'nomessage 'noerror
;; NB: loading the ess-module slows down emacs' load time a lot, even if ESS is compiled.

;; ##########################################################################



(when window-system
  (mwheel-install) ;; enable wheelmouse
  (setq focus-follows-mouse nil)
  (set-selection-coding-system 'compound-text-with-extensions) ;; fix X clipboard
  )

;; Set default modes for various file endings

(add-to-list 'auto-mode-alist '("Rakefile" . ruby-mode))
(add-to-list 'auto-mode-alist '("\\.rb\\'" . ruby-mode))

(add-to-list 'auto-mode-alist '("\\.mak\\'" . makefile-mode))
(add-to-list 'auto-mode-alist '("/\\Makefile[^/]*\\'" . makefile-mode)) ; Filename has "Makefile" somewhere in it
(add-to-list 'auto-mode-alist '("\\.txt\\'" . html-mode)) ; so that we get nice highlighting!
(add-to-list 'auto-mode-alist '("\\.tab\\'" . html-mode)) ; "text-mode" is so plain, so we use the nice HTML mode with highlighting
(add-to-list 'auto-mode-alist '("\\.lst\\'" . html-mode)) ; "text-mode" is so plain
(add-to-list 'auto-mode-alist '("\\.dat\\'" . html-mode)) ; "text-mode" is so plain


;;(defun insert-literal-char (arg)
;;  "Note: Ctrl-Q inserts literal *typeable* characters."
;;  "Insert a character into a buffer by specifying its ascii code"
;;  (interactive "nEnter decimal value of ASCII chracter to insert: ") ;; note: 9 is tab (so if arg=9, it inserts a tab)
;;  (insert (format "%c" arg))
;;)

;; Turn off all the annoying distractions
(if (fboundp 'tool-bar-mode) (tool-bar-mode -1))
(if (fboundp 'menu-bar-mode) (menu-bar-mode -1))
;; - http://www.cabochon.com/~stevey/blog-rants/effective-emacs.html


(show-paren-mode          t) ;; highlight matching parentheses
(auto-compression-mode    1) ;; transparently reads  .gz files
(global-font-lock-mode    t) ;; turn on syntax highlighting

(delete-selection-mode 1) ; delete seleted text when typing
;(global-linum-mode 1) ; always show line numbers
;(global-visual-line-mode 1) ; wrap lines at word boundaries instead of by characters

(setq
 inhibit-startup-message   t
 frame-title-format        (list "Emacs: %f") ;; Set titles: %f = filename, %b = buffer name
 icon-title-format         "Emacs - %b"	;; see above line
 truncate-partial-width-windows  nil ; <-- "nil" allows word wrap to function  with vertically-split (Ctrl-X 3) windows. Use "toggle-truncate-lines" to change it
 scroll-step               1 ; Scroll only this many lines when we hit the bottom
 scroll-conservatively     1
 show-trailing-whitespace  t
 tab-width                 8 ;; <-- tab indent width!
 indent-tabs-mode          t
 default-fill-column       80 ;; <-- when you "meta-q" to fit text, what line widths are used for wrapping?
 transient-mark-mode       t ;; show when we mark text for copying/selection, BUT also makes the selection disappear after one use. (if unset, we can use ctrl-space ctrl-space)

 search-highlight          t
 isearch-lazy-highlight    t ;; t = "show all search results"
 query-replace-highlight   t ;; <-- highlight during searching
 abbrev-case-fold-search   t

 blink-matching-paren      t
 blink-matching-paren-distance 999000
 line-number-mode          t ; Show current line number in status bar
 column-number-mode        t ; Show current column number in status bar
 which-func-mode           t
 kill-whole-line           nil ; Removes newline with Ctrl-K, if cursor is at beginning
 visible-bell              nil ; nil = No screen flashing
 comint-input-ring-size    1000 ; default "history" is just 32
 default-major-mode       'text-mode ; default mode
 require-final-newline     t ; always end a file with a newline
 inhibit-default-init      t ; Don't auto-load additional config files
 next-line-add-newlines    nil ; don't insert newlines when the cursor moves beyond the file's end
 sentence-end-double-space nil
 colon-double-space        nil
 sentence-end              "[.?!][]\"')]*\\($\\|\t\\| \\)[ \t\n]*"

 vc-follow-symlinks        t ; automatically follow symlinks to CVS-controlled files

 ;; delete-old-versions     t ;;  http://www.emacswiki.org/emacs/BackupDirectory
 kept-new-versions         6 ; BACKUPS
 kept-old-versions         2 ; BACKUPS
 backup-by-copying         t ; BACKUPS don't clobber symlinks (?)
 vc-handled-backend      nil
 vc-make-backup-files    nil ; BACKUPS
 make-backup-files         t ; BACKUPS
 auto-save-mode            t ; BACKUPS: Don't auto-save
 version-control         nil ; BACKUPS: use versioned backups

 cache-long-line-scans     t ; Makes emacs less slow when scrolling when there are very long lines
 show-paren-style  'expression ; highlight entire expression when a paren is selected
 ;; - http://users.tkk.fi/~rsaikkon/conf/dot.emacs
 )

(setq vc-handled-backends nil) ;; <-- don't load the incredibly slow "vc-hg" mercurial module

;; More about backups
(setq auto-save-mode t)
(setq backup-directory-alist '(("." . "~/.emacs-backups")))


(fset 'yes-or-no-p 'y-or-n-p) ; <-- Treats y/n as shortcuts for yes/no

;; <-- Note: Ctrl-Shift-O screws up the arrow keys for some reason!

(global-set-key [(meta o)] 'move-end-of-line)
(global-set-key [(meta u)] 'move-beginning-of-line)

(global-set-key (kbd "ESC <down>") 'scroll-up)
(global-set-key (kbd "ESC <up>") 'scroll-down)

(global-set-key [(meta a)] 'execute-extended-command)
;;(global-set-key "\M-o" 'forward-word)
;;(global-set-key "\M-u" 'backward-word)

(global-unset-key "\M-c")

(setq agw-defining-macro nil)
(global-set-key (kbd "M-\'")
		'(lambda () "Start keyboard macro" (interactive)
		   (if (equal agw-defining-macro nil)
		       (progn (setq agw-defining-macro (not agw-defining-macro))
			      (start-kbd-macro nil))
		     ;; else
		     (progn (setq agw-defining-macro (not agw-defining-macro))
			    (end-kbd-macro)))))

(global-set-key (kbd "M-;") '(lambda () "Run keyboard macro"   (interactive) (call-last-kbd-macro) (message "Repeated the last macro defined with M-[...")))



(global-set-key (kbd "M-1") '(lambda () "Close other windows"   (interactive) (delete-other-windows) (message "delete-other-windows: Closed all windows besides the active one.")))

(global-set-key (kbd "M-2") '(lambda () "Split window vertically"   (interactive) (split-window-vertically) (message "split-window-vertically: Split the window into top/bottom panes.")))

(global-set-key (kbd "M-3") '(lambda () "Split window horizontally"   (interactive) (split-window-horizontally) (message "split-window-horizontally: Split the window into left/right panes.")))

(global-set-key (kbd "M-%") 'set-variable)
(global-set-key (kbd "M-5") 'describe-variable)
(global-set-key (kbd "M-6") 'describe-key)

(global-set-key (kbd "M-7") '(lambda () "Delete Trailing Whitespace"   (interactive) (delete-trailing-whitespace) (message "delete-trailing-whitespace: Deleted trailing whitespace (if any)")))

(global-set-key (kbd "M-8") '(lambda () "Adds a shell to the bottom quarter of this window."
			       (interactive) (split-window-vertically) (other-window 1) (split-window-vertically) (delete-window) (shell) (rename-buffer "Shell-primary-buffer") (other-window 1) (message "Set up this window with a shell at the bottom.")))

(global-set-key (kbd "M-9") 'previous-buffer)
(global-set-key (kbd "M-0") 'next-buffer)

(global-set-key (kbd "M--") '(lambda () "Close this window"   (interactive) (delete-window) (message "Just ran delete-window (M-0)! (Note: the buffer is still active. Use <C-x b> to find it)")))


(global-set-key [(meta m)] '(lambda () "Next pane..." (interactive) (other-window 1)))

(global-set-key [(shift meta m)] '(lambda () "Previous pane..." (interactive) (other-window -1)))

;;(global-unset-key "\C-v")
;;(global-unset-key "\M-v") ;; now it's meta-u and meta-o


(global-unset-key [(control x) (m)])

(global-unset-key [insert])	  ; Disable the "insert" key
(global-set-key [insert] (function (lambda () (interactive) (message "The insert key was DISABLED in the ~/.emacs file."))))

;; only set the html-mode-map if it isn't ALREADY defined
(when (not (fboundp 'html-mode-map))
  (setq html-mode-map (make-sparse-keymap))
  (define-key html-mode-map "\t" 'tab-to-tab-stop)
  (define-key html-mode-map (kbd "TAB") 'self-insert-command) ; tab inserts ONE TAB ONLY AND DOES NOT INDENT in this mode
  )

(define-key text-mode-map (kbd "TAB") 'self-insert-command) ; tab inserts ONE TAB ONLY AND DOES NOT INDENT in text mode

(define-key text-mode-map (kbd "M-s") 'save-buffer) ; tab inserts ONE TAB ONLY AND DOES NOT INDENT in text mode


;; FIX THE BACKSPACE AND DELETE KEYS
;; (finally they work like they're supposed to)
;;(global-set-key '(meta [return]) '[return])
(normal-erase-is-backspace-mode 0)
(global-set-key [delete] 'delete-char)
(global-set-key [kp-delete] 'delete-char)
(global-set-key [backspace] 'backward-delete-char-untabify)
(global-set-key [C-backspace] 'backward-delete-char-untabify)
(global-set-key [S-backspace] 'kill-word)
(define-key isearch-mode-map [backspace] 'isearch-delete-char) ; lets you backspace in a search string (ctrl-S now works right)
;;(define-key isearch-mode-map "\C-h" 'isearch-delete-char) ; lets you delete in a search string (ctrl-S now works right)
;;(global-set-key "\C-h" 'backward-delete-char-untabify)
;;(global-set-key "\M-\C-h" 'backward-delete-char)

(global-set-key "\M-s" 'save-buffer)
(global-set-key "\M-z" 'undo)
(global-set-key "\M-e" 'recentf-open-files)
(recentf-mode 1)

;;         ^  Navigation Keys (ijkl)
;;         |  (Inverted-T navigation)
;;         I
;;   <- J  K  L -->
;;         |
;;         v
(global-set-key "\M-i" 'previous-line)
(global-set-key "\M-k" 'next-line)
(global-set-key "\M-j" 'backward-char)
(global-set-key "\M-l" 'forward-char)

(global-set-key "\M-J" '(lambda () "Moves left quickly"   (interactive) (MoveCursorWithoutWrapping -8)))
(global-set-key "\M-L" '(lambda () "Moves right quickly"    (interactive) (MoveCursorWithoutWrapping 8)))
(global-set-key "\M-I" '(lambda () "Moves the cursor up quickly"  (interactive) (previous-line 5)))
(global-set-key "\M-K" '(lambda () "Moves the cursor down quickly" (interactive) (next-line 5)))

(global-set-key (kbd "M-/") 'comment-dwim) ;; or-uncomment-region)

(global-set-key "\M-\t" 'dabbrev-expand)


;; Tab-related commands
;;(global-set-key "\t" 'indent-region)
;;(global-set-key [(shift \t)] '(insert (format "%c" "\t")))


(global-set-key (kbd "M-.") 'dabbrev-expand)
(global-set-key (kbd "M-SPC") 'self-insert-command) ;; Causes meta-space to also make a space!

(global-set-key (kbd "M-!") 'hs-hide-block) ;; Code folding!
(global-set-key (kbd "M-@") 'hs-show-block)

(setq agw-show-all t) ;; <-- default: show all code-folding folds
(global-set-key (kbd "M-#") '(lambda () "Hidw/show all code folds."
			       (interactive)
			       (if (equal agw-show-all nil)
				   (progn (hs-show-all) (message "Show all code folds"))
				 (progn (hs-hide-all) (message "Hide all folds")))
			       (setq agw-show-all (not agw-show-all))))



(global-unset-key (kbd "M-g"))
(global-set-key (kbd "M-G") '(lambda () "Mercurial commit from directly within emacs" 
			       (interactive) 
			       (shell-command
				(concat 
				 "hg commit "
				 " --user " "'" (getenv "USER") "'"
				 " --message \"Commit from within emacs\""))))

;;(global-set-key (kbd "M-G") '(lambda () "Git commit" (interactive) (shell-command "git commit -a -m \"Commit from within emacs\"")))


;;;  Bind my-command to f1
;;(global-set-key 'f1 'my-command)

;;;  Bind my-command to Shift-f1
;;(global-set-key '(shift f1) 'my-command)

;;; Bind my-command to C-c Shift-f1
;;(global-set-key '[(control c) (shift f1)] 'my-command)

;;; Bind my-command to the middle mouse button.
;;(global-set-key 'button2 'my-command)

;;(global-set-key "\M-O" 'end-of-buffer) ; <-- for some reason, this causes the arrow keys to stop working
(global-set-key "\C-t" 'toggle-truncate-lines)
(global-set-key "\C-o" 'find-file)
(global-set-key "\C-j" 'goto-line)
(global-set-key "\M-r" 'query-replace) ; less annoying than meta-shift-5
(global-set-key [(control meta r)] 'query-replace-regexp) ; less annoying than meta-shift-5
(global-set-key "\C-s" 'isearch-forward)
(global-set-key "\C-\\" 'indent-region)
(global-set-key "\M-\C-\\" 'indent-region)


;; Redo syntax coloring when it gets screwed up
(global-set-key (kbd "C-4")
		'(lambda () "Redoes the syntax coloring..." (interactive)
		   (global-font-lock-mode nil)
		   (global-font-lock-mode t)
		   (message "Redid the syntax coloring (Ctrl-m)")
		   ))


(defun ChangeTabWidthBy(AMT)
  "AGW: Change tab width by some amount."
  (let ((newWidth (+ tab-width AMT)))
    (set-variable 'tab-width (cond
			      ((> newWidth 2) newWidth)
			      (t 2)))
    (message (concat "Changed the tab width to "
		     (int-to-string tab-width)))))

(global-set-key "\M-T" '(lambda () "Reduce tab width"   (interactive) (ChangeTabWidthBy -2)))
(global-set-key "\M-t" '(lambda () "Increase tab width"   (interactive) (ChangeTabWidthBy 2)))


(defun MoveCursorWithoutWrapping(MOVE_AMT)
  "AGW: Move the cursor horizontally, but don't wrap past the end of a line"
  (interactive)
  (let ((MOVE_TO_COL (+ MOVE_AMT (current-column))))
    (cond ((< MOVE_TO_COL 0) (move-to-left-margin))
	  (t (move-to-column MOVE_TO_COL)))))

;; COLORS!

(defun agwColor (type fore back)
  (if (not (null fore)) (set-face-foreground type fore))
  (if (not (null back)) (set-face-background type back)))
;;(copy-face 'bold 'font-lock-function-name-face)

;; ======== <HIGHLIGHT THE CURRENT LINE>: ====
;;(global-hl-line-mode 1) ; line highlight / highlight line with the active cursor on it
;;(set-face-background 'hl-line "gray10")
;; ======== </HIGHLIGHT THE CURRENT LINE> ====

(set-cursor-color CURSOR_BG)
(agwColor 'cursor CURSOR_TEXT CURSOR_BG)
;;(set-background-color BACKGROUND) ;; Set emacs background color
;;(set-foreground-color FOREGROUND) ;; Set emacs background color

;; MODE LINE COLORS
(agwColor 'modeline ACTIVE_MODELINE ACTIVE_MODELINE_TEXT)
(agwColor 'mode-line-inactive INACTIVE_MODELINE INACTIVE_MODELINE_TEXT)
(agwColor 'mode-line-highlight "green" "red")
(agwColor 'mode-line-buffer-id BUFFER_FILENAME_TEXT BUFFER_FILENAME_BG)
(agwColor 'minibuffer-prompt MINIBUFFER_TEXT MINIBUFFER_BG)


;; SEARCH HIGHLIGHTING
(agwColor 'isearch SEARCH_RESULT_HIGHLIGHT_TEXT SEARCH_RESULT_HIGHLIGHT_BG)
(set-default 'query-replace 'isearch)
(agwColor 'lazy-highlight SEARCH_RESULT_OTHER_TEXT SEARCH_RESULT_OTHER_BG)


(agwColor 'completions-common-part "cyan" "purple")

(agwColor 'show-paren-match "yellow" "red")
(agwColor 'show-paren-mismatch ERR_TEXT ERR_BG)

(agwColor 'trailing-whitespace "black" "magenta")

(agwColor 'nobreak-space ERR_TEXT ERR_BG)
(agwColor 'escape-glyph ERR_TEXT ERR_BG) ; Highlights the '\' or '^' escape character


(custom-set-faces
  ;; custom-set-faces was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(font-lock-comment-face ((((class color)) (:bold t :underline nil :foreground "blue" :background nil)))))
;(agwColor 'font-lock-comment-face COMMENT_FOREGROUND COMMENT_BACKGROUND)
(set-default 'font-lock-comment-delimiter-face 'font-lock-comment-face)

(agwColor 'region SELECTED_TEXT SELECTED_BG) ; For a selected region

;; Programming Syntax Styles
(agwColor 'font-lock-function-name-face FUNCTION_TEXT FUNCTION_BG) ; For function names
(agwColor 'font-lock-keyword-face KEYWORD_TEXT KEYWORD_BG)

(agwColor 'font-lock-type-face TYPE_TEXT TYPE_BG) ; Face name to use for type and class names
(agwColor 'font-lock-string-face STRING_TEXT STRING_BG)

(agwColor 'font-lock-variable-name-face VARIABLE_TEXT VARIABLE_BG) ; Font Lock mode face used to highlight variable names.
(agwColor 'font-lock-constant-face "orange" VARIABLE_BG)

(agwColor 'font-lock-negation-char-face "red" "white")

(agwColor 'font-lock-preprocessor-face FUNCTION_TEXT FUNCTION_BG)

(agwColor 'font-lock-regexp-grouping-backslash "black" "blue")
(agwColor 'font-lock-regexp-grouping-construct "black" "blue")

(agwColor 'font-lock-warning-face ERR_TEXT ERR_BG)

(agwColor 'font-lock-builtin-face BUILTIN_TEXT BUILTIN_BG) ; For built-in functions

(agwColor 'font-lock-doc-face DOC_TEXT DOC_BG) ; Face name to use for documentation.

;(add-hook 'font-lock-mode-hook 'highlight-hard-spaces)
;; (add-hook 'font-lock-mode-hook 'highlight-tabs) ;; <-- uncomment if you want tabs ALWAYS highlighted
;font-lock-mode-hook
;(add-hook 'makefile-mode-hook 'highlight-trailing-whitespace) ;; <-- show trailing whitespaces in Makefile mode


(add-hook 'c-mode-common-hook   'hs-minor-mode)
(add-hook 'emacs-lisp-mode-hook 'hs-minor-mode)
(add-hook 'java-mode-hook       'hs-minor-mode)
(add-hook 'lisp-mode-hook       'hs-minor-mode)
(add-hook 'perl-mode-hook       'hs-minor-mode)
(add-hook 'sh-mode-hook         'hs-minor-mode)
(add-hook 'ess-mode-hook        'hs-minor-mode)

(add-hook 'shell-mode-hook 'ansi-color-for-comint-mode-on)


(if should-load-ess
    (progn
      (add-hook 'ess-mode-hook  ;; R-mode-hook r-mode-hook r mode <-- should be ess-mode-hook
		(progn (define-key ess-mode-map "\M-\t" 'dabbrev-expand)  ;; Make meta-tab do the normal expansion even in ESS mode
		       (define-key ess-mode-map "_" nil)          ;; no smart underscores!
		       (define-key inferior-ess-mode-map "_" nil) ;; no smart underscores!
		       ;;		 (define-key ess-mode-map "\t" 'self-insert-command)
		       ))
      
      ;; note: lambda and progn are different somehow!!!
      (add-hook 'ess-mode-hook
		(lambda ()
		  (ess-set-style 'BSD)
		  (setq ess-indent-level 5)
		  (setq ess-fancy-comments 'nil)
	      ))
      ))

(add-hook 'linum-before-numbering-hook
	  (lambda ()
	    (let ((numDigitsInTotalNumLines (length (number-to-string (count-lines (point-min) (point-max)))))) 
	      (setq linum-format (concat "%" (number-to-string numDigitsInTotalNumLines) "d "))))) ;; makes a result like "%5d " or "%4d " for formatting, which then is passed to the command (setq linum-format "%5d ")))

(add-hook 'makefile-mode-hook
	  (lambda ()
	    (define-key makefile-mode-map "\M-\t" 'dabbrev-expand)  ;; Make meta-tab do the normal expansion even in Make mode
	    ))

;(add-hook 'makefile-mode-hook 'highlight-tabs)

(add-hook 'text-mode-hook
	  (lambda ()
	    (define-key text-mode-map "\M-\t" 'dabbrev-expand)  ;; Make meta-tab do the normal expansion even in text mode (default is the terrible ispell mode)
	    ))


(font-lock-add-keywords
 'perl-mode
 '(
   ("\\<\\(DEBUG\\|Debug\\|debug\\|DEBUGGING\\|Debugging\\|debugging\\)\\>" 1 font-lock-warning-face prepend)
   ("-" . font-lock-builtin-face)
   )
 )

(font-lock-add-keywords
 'makefile-mode
 '(
   ("\\<\\([0-9]+.tmp\\)" 1 agwTempFilenameFace prepend)
   ("\\<\\(DEBUG\\|Debug\\|debug\\|DEBUGGING\\|Debugging\\|debugging\\)\\>" 1 font-lock-warning-face prepend)
   ("\\<\\([-0-9A-Za-z_]+.tmp\\)" 1 agwTempFilenameFace keep)
   ("\\<\\(arg_[-0-9A-Za-z_]+\\)" 1 font-lock-builtin-face append)
   ("\\<\\(and\\|or\\|not\\)\\>"  1 font-lock-keyword-face keep)
   )
 )

(font-lock-add-keywords
 'python-mode
 '(
   ("\\<\\(pass\\|return\\|raise\\)\\>" 1 agwPythonPassFace keep)
   ("\\<\\(DEBUG\\|Debug\\|debug\\|DEBUGGING\\|Debugging\\|debugging\\)\\>" 1 font-lock-warning-face keep)
   )
 )


;;(defvar ALLOWED_R_VARIABLE_CHARS "[a-zA-Z0-9\\.]")

;;(defun testff () "test" (interactive)  (insert (concat "a" "b" ALLOWED_R_VARIABLE_CHARS)))




;; OVERRIDE and LAXMATCH are flags.  If OVERRIDE is t, existing fontification can be overwritten.  If `keep', only parts not already fontified are highlighted. If `prepend' or `append', existing fontification is merged with the new, in which the new or existing fontification, respectively, takes precedence. If LAXMATCH is non-nil, that means don't signal an error if there is no match for SUBEXP in MATCHER.

;; R mode special things (requires the ESS stuff to be loaded!)

					;(defvar zom "\"[abcdefghijV]\"")
					;(message "OK YOU SET THE VARIABLE TO:")
					;(message zom)
					;(message "YEP THAT WAS IT")

;; Keyword highlight options: t (OVERRIDE) /  prepend (prioritize) / append / keep (keep existing coloring, if any)
(font-lock-add-keywords
 'ess-mode
 '(
   ;;   ("\\<\\(kv[a-zA-Z0-9\\.]*\\)\\($\\|[]-+~` 	<>=,;:(){}%*!@#$^&\\/\'\"]\\)" 1 font-lock-constant-face append) ; anything that starts in kv
   ("\\<\\(stop\\|stopifnot\\|browser\\|options\\)\\>" 1 font-lock-warning-face keep)
   ("\\<\\(kv[a-zA-Z0-9\\.]*\\)\\($\\|[^a-zA-Z0-9\\.]\\)" 1 font-lock-constant-face keep) ; anything that starts in kv
   ("\\<\\(gv[a-zA-Z0-9\\.]*\\)\\($\\|[^a-zA-Z0-9\\.]\\)" 1 'agwMakeGlobalVarFace keep) ; anything that starts in gv

   ;;("^\\([}].*\\)"      1 'agwIndent1Face t) ; function-ending (line-starting) brace
   
   ("^\\([ ]\\{1,1\\}\\)"      1 'agwIndent1Face t) ; line-starting tab
   ("^    \\([ ]\\{1,1\\}\\)" 1 'agwIndent2Face t) ; line-starting tab
   ("^         \\([ ]\\{1,1\\}\\)" 1 'agwIndent3Face t) ; line-starting tab
   ("^              \\([ ]\\{1,1\\}\\)" 1 'agwIndent4Face t) ; line-starting tab
   ("^                   \\([ ]\\{1,1\\}\\)" 1 'agwIndent3Face t) ; line-starting tab

   ("\\(\<\<-\\)" 1 'agwMakeGlobalVarFace keep) ; the <<- global assignment operator
   ("\\(\<-\\)" 1 'agwPositiveNumberFace keep) ; the <- regular assignment operator
   ("\\<\\(assert\.agw.*\\)" 1 'agwAssertionFace t) ; anything that ends in Vec
   ("\\<\\([a-zA-Z0-9_\\.]*Vec\\)\\($\\|[^a-zA-Z0-9\\.]\\)" 1 'agwArrayFace keep) ; anything that ends in Vec
   ("\\<\\([a-zA-Z0-9_\\.]*List\\)\\($\\|[][-+~` 	<>=,;:(){}%*!@#$^&\\/\'\"]\\)" 1 'agwListFace keep) ; anything that ends in List
   ("\\<\\([a-zA-Z0-9_\\.]*Hash\\)\\($\\|[]-+~` 	<>=,;:(){}%*!@#$^&\\/\'\"]\\)" 1 'agwListFace keep) ; anything that ends in Hash
   ("\\<\\([a-zA-Z0-9_\\.]*Mat\\)\\($\\|[]-+~` 	<>=,;:(){}%*!@#$^&\\/\'\"]\\)" 1 'agwListFace keep) ; anything that ends in Hash
   ("\\<\\([a-zA-Z0-9_\\.]*\\.vec\\)\\($\\|[^a-zA-Z0-9\\.]\\)" 1 'agwArrayFace keep) ; anything that ends in Vec
   ("\\<\\([a-zA-Z0-9_\\.]*\\.list\\)\\($\\|[][-+~` 	<>=,;:(){}%*!@#$^&\\/\'\"]\\)" 1 'agwListFace keep) ; anything that ends in List
   ("\\<\\([a-zA-Z0-9_\\.]*\\.hash\\)\\($\\|[]-+~` 	<>=,;:(){}%*!@#$^&\\/\'\"]\\)" 1 'agwListFace keep) ; anything that ends in Hash
   ("\\<\\([a-zA-Z0-9_\\.]*\\.mat\\)\\($\\|[]-+~` 	<>=,;:(){}%*!@#$^&\\/\'\"]\\)" 1 'agwListFace keep) ; anything that ends in Hash
   ("\\([=]=======.*\\)" 1 'agwCustomDoubleLineFace t) ;; <-- eight '=' in a row means "highlight this line in a visually obvious manner"
   )
 )

(add-hook
 'after-change-major-mode-hook
 '(lambda ()
    (font-lock-add-keywords
     nil
     '(
;;       ("\\([0-9]+\\)" 1 'agwLineNumberFace t) ; the <- regular assignment operator
 ;      ("\\<\\([+]?[0-9\\.]+\\(?:\\.[0-9]+\\)?\\(?:[eE][-+]?[0-9]+\\)?\\)\\>" 1 'agwPositiveNumberFace keep)
       ("\\(-[0-9\\.]+\\(?:\\.[0-9]+\\)?\\(?:[eE][-+]?[0-9]+\\)?\\)" 1 'agwNegativeNumberFace keep)
       ("\\<\\(NaN\\|NA\\|ND+\\)\\>" 1 font-lock-warning-face keep)
       ("\\(.*[=]=======.*\\)" 1 'agwCustomDoubleLineFace t) ;; <-- eight '=' in a row means "highlight this line in a visually obvious manner"
       ))))

;; FUNCTIONS
(defun dos-to-unix-convert-line-endings () (interactive) (goto-char (point-min))
  (while (search-forward "\r" nil t) (replace-match "")))
(defun mac-to-unix-convert-line-endings () (interactive) (goto-char (point-min))
  (while (search-forward "\r" nil t) (replace-match "\n")))

;;make-comment-invisible
(defun hide-comments ()
  (interactive) (custom-set-faces '(font-lock-comment-face ((((class color)) (:foreground "black" :background "black")))))
  (message "Hide comments: Made comments invisible (set their color to the background color.")
  )

;;make-comment-visible
(defun show-comments ()
  (interactive) (custom-set-faces '(font-lock-comment-face ((((class color)) (:foreground "dark red" :background "black")))))
  (message "Show comments: Made comments visible again.")
  )

;;(global-set-key "\M-9" 'hide-comments)
;;(global-set-key "\M-(" 'show-comments)

;;(global-set-key "\C-x w" 'kill-current-buffer-and-window)


;;(setq scroll-preserve-screen-position 1)

(setq make-backup-files  nil) ;; Don't make those #scratch-042# backup files!
(auto-save-mode          nil) ;; Don't auto save!!!!!

;; ================= BACKUP FILES AND AUTOSAVING =================

(iswitchb-mode t) ;; enhances buffer switching (C-x b)

;;(message "My .emacs loaded in %ds" (destructuring-bind (hi lo ms) (current-time) (- (+ hi lo) (+ (first *emacs-load-start*) (second *emacs-load-start*)))))


;;(defun CustomMoveByWord(X)
;;  "AGW: Move the cursor by X words horizontally."
;;  (interactive)
;;  (skip-chars-forward "a-zA-Z "))

;;  (forward-word X))

;;(global-set-key "\M-p" 'delete-trailing-whitespace)

;;(set-comment-column 15)

;;(define-key minibuffer-local-map (kbd "C-<tab>") 'dabbrev-expand)

;; Hide lines that don't match a given search string
(setq has-hidelines nil)
;(if (system-type-is-darwin)
;    (progn (autoload 'hide-lines "hide-lines" "Hide lines based on a regexp" t) ;; ~/.emacs.d/hide-lines.el
;	   (require 'hide-lines nil t)
;	   (require 'hidesearch nil t)
;	   (setq has-hidelines t)))


(global-set-key [(control meta l)] '(lambda () (interactive) (hidesearch) (message "Only showing matching lines: use Ctrl-G to show all lines again."))) ;; control g now shows the invisibles, as well as generally cancelling


(global-set-key [(control g)] '(lambda () "Cancels normal operations, and also cancels out of search-and-hide-lines mode."
				 (interactive)
				 (progn (if has-hidelines (show-all-invisible))
					(message "Super cancel!")
					(keyboard-quit)
					)))

(global-set-key [(shift meta f)] 'grep-find)
(global-set-key [(meta control s)] 'grep-find)

;(global-set-key [(meta c)] 'kill-ring-save)
;(global-set-key [(meta x)] 'kill-region)
;(global-set-key [(meta v)] 'yank)
;(global-set-key [(shift meta v)] 'yank-pop)


;;(split-window-horizontally)   ;; want two windows at startup
;;(other-window 1)              ;; move to other window
;;(shell)                       ;; start a shell
;;(rename-buffer "shell-first") ;; rename it
;;(other-window 1)              ;; move back to first window

(global-set-key (kbd "RET") 'newline) ;; Make sure RET inserts a newline!


(defun acount(string)
  "Counts the number of occurences of STRING in the buffer."
  (interactive (list (read-string "Count: ")))
  (let ((point (point))
	(counter 0)
	(myStr   nil)
	(mat nil))
    (beginning-of-buffer)
    (while (re-search-forward string (point-max) t)
      (setq ps (match-beginning 0))
      (setq pe (match-end 0))
      (setq counter (+ 1 counter))
      (setq myStr (concat myStr (match-string 0) "\n"))
      ;(setq myStr (concat myStr "<" (thing-at-point 'symbol) ">=" (format "<%d," ps) (format " %d>  " pe) (match-string 0) "   :: " ))
      )
    (goto-char point)
    (message (format "%d occurences" counter))
    (message myStr)
    ))

(defun am()
  "Fun!"
  (interactive)
  (acount "^[^# 	]+.*\\(=\\|<-\\)[ ]*function[ ]*(")
  ;(acount "^.*\\(=\\|<-\\)[ ]*function[ ]*(")
  ) ;(count-occurences "alias"))


(defun reload-config () 
  "Runs load-file on ~/.emacs" 
  (interactive)
  (load-file "~/.emacs")
  (message "Loaded the config file again..."))

(defun buffer-exists (bufname)   (not (eq nil (get-buffer bufname))))

;; Close the loathsome scratch buffer!
(if (buffer-exists "*scratch*")  (kill-buffer "*scratch*"))
(if (buffer-exists "*Messages*")  (kill-buffer "*Messages*"))

;; ##

(custom-set-variables
  ;; custom-set-variables was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(case-fold-search t)
 '(read-buffer-completion-ignore-case t)
 '(read-file-name-completion-ignore-case t)
 '(completion-ignore-case t)
 '(apropos-do-all t))


;; (when (require 'browse-kill-ring nil t)

;;     ;; string separating entries in the `separated' style
;;     (setq browse-kill-ring-separator
;;           "\n--separator------------------------------")

;;     ;; temporarily highlight the inserted `kill-ring' entry
;;     (setq browse-kill-ring-highlight-inserted-item t)

;;     ;; face in which to highlight the `browse-kill-ring-separator'
;;     (defface separator-face '((t (:foreground "Blueviolet" :weight bold))) nil)
;;                                         ; slate gray
;;     (setq browse-kill-ring-separator-face 'separator-face)

;;     ;; use `M-y' to invoke `browse-kill-ring'
;;     (browse-kill-ring-default-keybindings))


