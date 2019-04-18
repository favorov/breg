#!/bin/bash
CHAINS=(IGH IGL IGK)
SHARES=(0 0.05)
for chain in ${CHAINS[@]}; do
	for share in ${SHARES[@]}; do
		echo ./combine_clones --chain ${chain} --mismatches-share ${share} --clones-folder ../clones/
		./combine_clones --chain ${chain} --mismatches-share ${share} --clones-folder ../clones/
		echo done...
	done
done
