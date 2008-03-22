#!/usr/bin/env ruby

require 'tempfile'
require 'optparse'

centroid_fold="centroid_fold"
simfold_pf="simfold_pf"
th=1e-5
opts=""
params=""
OptionParser.new do |opt|
  opt.on('-g', '--gamma gamma', 'weight of base-pairs') {|v| opts+=" -g #{v}"}
  opt.on('--mea', 'MEA estimators') {|v| opts+=" --mea" }
  opt.on('-t', '--threshold th', 'threshold of posteriors') {|v| th=v.to_f}
  opt.on('--centroid_fold path', 'exec path of centroid_fold') {|v| centroid_fold = v}
  opt.on('--simfold_pf path', 'exec path of simfold_pf') {|v| simfold_pf = v}
  opt.on('--params params',
         'specify a parameter file for simfold') {|v| params = "-p #{v}"}
  opt.parse!(ARGV)
end

def load_aln(f)
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
  desc.map{|d| seq[d].delete('-')}
end

ARGV.each do |f|
  seqs=load_aln(f)
  p = seqs.map do |seq|
    seq_c=seq.split('')
    bp=[]
    seq_c.each_with_index{|c,i| bp.push(["#{i+1} #{c}"])}
    IO.popen("#{simfold_pf} -s #{seq} #{params} -t #{th}", "r") do |io|
      cur=-1
      while o=io.gets
        if o=~/^\d/
          left,right,prob=o.chomp.split
          left=left.to_i
          right=right.to_i
          bp[left].push("#{right+1}:#{prob}")
        end
      end
    end
    fh=Tempfile.open(File.basename(f))
    bp.each{|e| fh.puts e.join(" ")}
    fh.close
    fh
  end
  pos=p.map{|q| q.path}.join(' ')
  system("#{centroid_fold} #{opts} --aux #{f} #{pos}")
end
