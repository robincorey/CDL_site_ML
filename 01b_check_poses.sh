get_dist () {
out=$1/$2/$3
if [[ ! -f $out/site_CL.eq.xvg ]]
then
        gmx distance -f $out/eq.xtc -s $out/eq.tpr -oxyz $out/site_CL.eq.xvg -select 'cog of group "Site" plus cog of group "CARD_HG"' -rmpbc no -pbc no -n $out/site.ndx >& $out/out_files/out_dist_eq
fi
read -r t0 x0 y0 z0 <<<$(grep " 0.000" $out/site_CL.eq.xvg)
read -r t50 x50 y50 z50 <<<$(grep " 50000" $out/site_CL.eq.xvg)
read -r t100 x100 y100 z100 <<<$(grep " 100000" $out/site_CL.eq.xvg)
echo $1 $2 $3 `echo "scale = 4; sqrt($x0^2 + $y0^2)" | bc` `echo "scale = 4; sqrt($x50^2 + $y50^2)" | bc` `echo "scale = 4; sqrt($x100^2 + $y100^2)" | bc`
}

echo pdb site pose t0 t50 t100

for site in `ls 1*/*/*/eq.xtc`
do
        get_dist `echo $site | cut -f1 -d'/'` `echo $site | cut -f2 -d'/'` `echo $site | cut -f3 -d'/'`
done
