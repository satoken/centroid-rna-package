#!/usr/bin/env ruby

require 'tempfile'
require 'optparse'

centroid_fold="centroid_fold"
contrafold="contrafold"
th=1e-5
opts=""
viterbi=false
params=""
OptionParser.new do |opt|
  opt.on('-g', '--gamma gamma', 'weight of base-pairs') {|v| opts+=" -g #{v}"}
  opt.on('--mea', 'MEA estimators') {|v| opts+=" --mea" }
  opt.on('--viterbi', 'viterbi estimators') {|v| viterbi=true }
  opt.on('-t', '--threshold th', 'threshold of posteriors') {|v| th=v.to_f}
  opt.on('--centroid_fold path', 'exec path of centroid_fold') {|v| centroid_fold = v}
  opt.on('--contrafold path', 'exec path of contrafold') {|v| contrafold = v}
  opt.on('--params params',
         'specify a parameter file for contrafold') {|v| params = "--params #{v}"}
  opt.parse!(ARGV)
end

ARGV.each do |f|
  open(f) do |fh|
    fh.gets('>')
    while l=fh.gets('>')
      fa=Tempfile.open(File.basename(f))
      l=l.chomp('>').split("\n")
      fa.puts '>'+l.shift
      l.each {|ll| fa.puts ll unless ll=~/[\(\)]/ }
      fa.close
      if viterbi
        system("#{contrafold} predict --viterbi #{fa.path}")
      else
        p=Tempfile.open(File.basename(f))
        p.close
        system("#{contrafold} predict #{params} --posteriors #{th} #{p.path} #{fa.path}")
        system("#{centroid_fold} #{opts} --aux #{fa.path} #{p.path}")
      end
    end
  end
end
