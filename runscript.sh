for i in {1..3}
do
  src/lda -est -alpha 0.05 -beta 0.01 -ntopics 20 -niters 1000 -twords 50 -fname kick-$i-cons -confile ../../files/kickstarterconstraints -dfile ../../files/kickstarter/headkickstarter.txt -savestep 1000 > kick-$i-cons.txt
done
