const express = require("express");
const { exec } = require("child_process");
const fs = require("fs");
// import code from "./code";
const app = express();
const PORT = 3000;

app.set("view engine", "ejs");
app.use(express.static("public"));

app.use(express.urlencoded({ extended: true }));

// Parse JSON bodies (as sent by API clients)
app.use(express.json());

app.get(`/`, (req, res) => {
  res.redirect('/work');
});

app.post("/getCode", (req, res) => {
  // console.log(req.body);
  // res.send(req.body);
  const input = `${req.body.n_value} ${req.body.x_value} ${req.body.y_value} ${req.body.e_value}`;
  fs.writeFileSync("input.txt", input);
  exec("g++ algo_code.cpp -o algo_code.exe", (error, stdout, stderr) => {
    if (error) {
      console.error(`Compilation error: ${error.message}`);
      res.send("error");
      return;
    }

    console.log("C++ code compiled successfully.");

    // Run C++ executable
    exec("algo_code", (error, stdout, stderr) => {
      if (error) {
        console.error(`Execution error: ${error.message}`);
        return;
      }
      const filePath = "output.txt"; // Adjust the path to your text file
      const content = fs.readFileSync(filePath, "utf-8");

      // Render the EJS template and pass the content to it
      res.render("result.ejs", { content });
      // res.render("result.ejs", {content});
      // res.send("done");
      // res.send(stdout);
      // console.log(`C++ program output:\n${stdout}`);
    });
  });
});

app.get(`/work`, (req, res) => {
  res.render("work.ejs");
});

app.listen(PORT, () => {
  console.log(`Listening at port ${PORT}`);
});
