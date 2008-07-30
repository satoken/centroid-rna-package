#!/usr/bin/env ruby

require 'tempfile'
require 'optparse'

$centroid_fold="centroid_fold"
$contrafold="contrafold"
$th=1e-5
$opts=""
$viterbi=false
$params=""
OptionParser.new do |opt|
  opt.on('-g', '--gamma gamma', 'weight of base-pairs') {|v| $opts+=" -g #{v}"}
  opt.on('--mea', 'MEA estimators') {|v| $opts+=" --mea" }
  opt.on('--viterbi', 'viterbi estimators') {|v| viterbi=true }
  opt.on('-t', '--threshold th', 'threshold of posteriors') {|v| $th=v.to_f}
  opt.on('--centroid_fold path', 'exec path of centroid_fold') {|v| $centroid_fold = v}
  opt.on('--contrafold path', 'exec path of contrafold') {|v| $contrafold = v}
  opt.on('--params params',
         'specify a parameter file for contrafold') {|v| $params = "--params #{v}"}
  opt.parse!(ARGV)
end

def predict_ss(f)
  open(f) do |fh|
    l=fh.gets('>')
    return false unless l=='>'
    while l=fh.gets('>')
      fa=Tempfile.open(File.basename(f))
      l=l.chomp('>').split("\n")
      fa.puts '>'+l.shift
      l.each {|ll| fa.puts ll unless ll=~/[\(\)]/ }
      fa.close
      if $viterbi
        system("#{$contrafold} predict --viterbi #{fa.path}")
      else
        p=Tempfile.open(File.basename(f))
        p.close
        system("#{$contrafold} predict #{$params} --posteriors #{$th} #{p.path} #{fa.path}")
        system("#{$centroid_fold} #{$opts} --aux #{fa.path} #{p.path}")
      end
    end
  end
  true
end

def aln2fa(f)
  desc=[]
  seq=Hash.new("")
  open(f) do |fh|
    fh.gets # skip header
    while l=fh.gets
      if l=~/^\S/
        n,s=l.chomp.split
        desc.push(n) unless seq.member?(n)
        seq[n] += s
      end
    end
  end

  return [] if desc.empty?
  
  tempfiles=[]
  desc.each do |n|
    fh=Tempfile.open(File.basename(f))
    fh.puts ">"+n
    fh.puts seq[n].delete('-')
    tempfiles.push(fh)
    fh.close
  end
  tempfiles
end

def predict_commom_ss(f)
  temp=aln2fa(f)
  return false if temp.empty?
  p = temp.map do |t|
    fh=Tempfile.open(File.basename(f))
    fh.close
    system("#{$contrafold} predict --posteriors #{$th} #{fh.path} #{t.path}")
    fh
  end
  pos=p.map{|q| q.path}.join(' ')
  system("#{$centroid_fold} #{$opts} --aux #{f} #{pos}")
  true
end

ARGV.each do |f|
  predict_ss(f) || predict_commom_ss(f)
end
