#! /bin/sh
pandoc $(printf '%s\n' *.md | sort) -o introduction_data.pdf