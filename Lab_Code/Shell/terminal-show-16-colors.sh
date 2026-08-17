#!/usr/bin/env bash
echo "Print the 'normal' and 'bold'-text 16-color terminal escape sequences in Bash:"
echo "Background | Foreground colors"
horizontalline="----------------------------------------------------------------------------"
echo $horizontalline
for((bg=40;bg<=47;bg++)); do
	for((bold=0;bold<=1;bold++)) do
		echo -en "\033[0m"" \\\033[${bg}m  | "
		if [ $bold == "0" ]; then echo -en "       " ;
		                     else echo -en "(BOLD) "
		fi
		for((fg=30;fg<=37;fg++)); do
			if [ $bold == "0" ]; then
				echo -en "\033[${bg}m\033[${fg}m [${fg}m  "
			else
				echo -en "\033[${bg}m\033[1;${fg}m [1;${fg}m"
			fi
		done
		echo -e "\033[0m"
	done
	echo $horizontalline
done

echo

