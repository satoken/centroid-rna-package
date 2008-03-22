#!/usr/bin/env ruby

require 'tempfile'
require 'optparse'

centroid_fold="centroid_fold"
simfold="simfold"
simfold_pf="simfold_pf"
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
  opt.on('--simfold path', 'exec path of simfold') {|v| simfold = v}
  opt.on('--simfold_pf path', 'exec path of simfold_pf') {|v| simfold_pf = v}
  opt.on('--params params',
         'specify a parameter file for simfold') {|v| params = "-p #{v}"}
  opt.parse!(ARGV)
end

ARGV.each do |f|
  open(f) do |fh|
    fh.gets('>')
    while l=fh.gets('>')
      l=l.chomp('>').split("\n")
      desc=l.shift
      seq=l.map{|ll| ll=~/[\(\)]/ ? '' : ll}.join('')
      if viterbi
        IO.popen("#{simfold} -s #{seq} #{params}", "r") do |io|
          while o=io.gets
            if o=~/MFE:\s+(\S+)\s+(-?[\d.]+)/
              puts ">"+desc
              puts seq
              puts "#{$1} (#{$2})"
            end
          end
        end
      else
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
        p=Tempfile.open(File.basename(f))
        bp.each{|e| p.puts e.join(" ")}
        p.close

        fa=Tempfile.open(File.basename(f))
        fa.puts ">"+desc
        fa.puts seq
        fa.close

        system("#{centroid_fold} #{opts} --aux #{fa.path} #{p.path}")
      end
    end
  end
end
