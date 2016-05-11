import java.io.IOException;
import java.util.StringTokenizer;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.fs.*;

import java.io.BufferedReader;
import java.io.InputStreamReader;

public class EditDistance {
  public static String X;
  public static String Y;
  public static int m, t;

  public EditDistance(String a, String b, int n) {
    X = a;
    Y = b;
    m = n;
  }

  public static void update(int x) {
    t = x;
  }

  public static class MyMapper
    extends Mapper<Object, Text, Text, LongWritable> {

    public Long sub(int i, int j) {
      return (X.charAt(i) == Y.charAt(j)) ? 0L : 1L;
    }

    public Long ins(long i, long j) {
      return (j * j) - (i * i);
    }

    public Long del(long i, long j) {
      return (j * j * j) - (i * i * i);
    }

    public void map(Object key, Text value, Context context
                   ) throws IOException, InterruptedException {
      String[] inputTokens = value.toString().split("\t");
      String[] inputKey = inputTokens[0].split(",");
      Text txt = new Text();
      txt.set(inputTokens[0]);

      long x = Long.parseLong(inputKey[0]);
      long y = Long.parseLong(inputKey[1]);
      Long val = Long.parseLong(inputTokens[1]);

      if (t % 2 == 1) {
        Text key1 = new Text();
        Text key2 = new Text();
        Text key3 = new Text();
        Long val1, val2, val3;

        if (x + 1 <= m && y + 1 <= m) {
          key1.set((x + 1) + "," + (y + 1));
          val1 = val + sub((int) x, (int) y);
          context.write(key1, new LongWritable(val1));
        }

        if (x + 1 <= m) {
          key2.set((x + 1) + "," + (y));
          val2 = val + ins(x, x + 1);
          context.write(key2, new LongWritable(val2));
        }

        if (y + 1 <= m) {
          key3.set((x) + "," + (y + 1));
          val3 = val + ins(y, y + 1);
          context.write(key3, new LongWritable(val3));
        }

        //context.write(txt, new LongWritable(val));
      } else {
        Text key1 = new Text();
        Long val1;

        if (x == y) {
          for (long i = x + 1; i <= m; i++) {
            key1.set((i) + "," + (y));
            val1 = val + ins(x, i);
            context.write(key1, new LongWritable(val1));
          }

          for (long j = y + 1; j <= m; j++) {
            key1.set((x) + "," + (j));
            val1 = val + ins(y, j);
            context.write(key1, new LongWritable(val1));
          }
        }
        context.write(txt, new LongWritable(val));
      }
    }
  }

  public static class MyReducer
    extends Reducer<Text, LongWritable, Text, LongWritable> {
    private LongWritable result = new LongWritable();

    public void reduce(Text key, Iterable<LongWritable> values,
                       Context context
                      ) throws IOException, InterruptedException {
      if (t % 2 == 0) {
        Long min = Long.MAX_VALUE;
        for (LongWritable val : values) {
          if (min > val.get())
            min = val.get();
        }
        result.set(min);
        context.write(key, result);
      } else {
        String[] inputKey = key.toString().split(",");
        long x = Long.parseLong(inputKey[0]);
        long y = Long.parseLong(inputKey[1]);
        int ind = t - (t / 2);

        if ((x == ind && y >= ind) || (y == ind && x >= ind)) {
          Long min = Long.MAX_VALUE;
          for (LongWritable val : values) {
            if (min > val.get())
              min = val.get();
          }
          result.set(min);
          context.write(key, result);
        }
      }
    }
  }

  public void run(String input, String outputFolder, int m) throws Exception {
    String output = outputFolder + System.nanoTime();
    for (int t = 1; t < 2 * m; t++) {
      update(t);
      Configuration conf = new Configuration();
      Job job = Job.getInstance(conf, "word count");

      job.setNumReduceTasks(40);
      job.setJarByClass(EditDistance.class);
      job.setMapperClass(MyMapper.class);
      job.setCombinerClass(MyReducer.class);
      job.setReducerClass(MyReducer.class);
      job.setOutputKeyClass(Text.class);
      job.setOutputValueClass(LongWritable.class);
      FileInputFormat.setInputDirRecursive(job, true);
      FileInputFormat.addInputPath(job, new Path(input));
      FileOutputFormat.setOutputPath(job, new Path(output));
      job.waitForCompletion(true);

      input = output + "/";
      output = outputFolder + System.nanoTime();
    }
  }
  public static void main(String[] args)  {
    try {
      Configuration conf = new Configuration();
      Path path = new Path(args[2]);
      FileSystem fs = FileSystem.get(path.toUri(), conf);
      FSDataInputStream inputStream = fs.open(path);
      BufferedReader br = new BufferedReader(new InputStreamReader(fs.open(path)));
      String x = br.readLine();
      
      path = new Path(args[3]);
      fs = FileSystem.get(path.toUri(), conf);
      inputStream = fs.open(path);
      br = new BufferedReader(new InputStreamReader(fs.open(path)));
      String y = br.readLine();
      fs.close();
      
      /*BufferedReader br = new BufferedReader(new FileReader("x.txt"));
      String x = br.readLine();
      br = new BufferedReader(new FileReader("y.txt"));
      String y = br.readLine();
      String x = "sleuth";
      String y = "health";*/
      
      int m = x.length();
      EditDistance ed = new EditDistance(x, y, m);
      ed.run(args[0], args[1], m);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}