#!/usr/bin/env ruby

require 'tempfile'
require 'optparse'

centroid_fold="centroid_fold"
contrafold="contrafold"
th=1e-5
opts=""
OptionParser.new do |opt|
  opt.on('-g', '--gamma gamma', 'weight of base-pairs') {|v| opts+=" -g #{v}"}
  opt.on('--mea', 'MEA estimators') {|v| opts+=" --mea" }
  opt.on('-t', '--threshold th', 'threshold of posteriors') {|v| th=v.to_f}
  opt.on('--centroid_fold path', 'exec path of centroid_fold') {|v| centroid_fold = v}
  opt.on('--contrafold path', 'exec path of contrafold') {|v| contrafold = v}
  opt.parse!(ARGV)
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

ARGV.each do |f|
  temp=aln2fa(f)
  p = temp.map do |t|
    fh=Tempfile.open(File.basename(f))
    fh.close
    system("#{contrafold} predict --posteriors #{th} #{fh.path} #{t.path}")
    fh
  end
  pos=p.map{|q| q.path}.join(' ')
  system("#{centroid_fold} #{opts} --aux #{f} #{pos}")
end
