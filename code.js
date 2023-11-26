const readline = require("readline");
const rl = readline.createInterface({
  input: process.stdin,
  output: process.stdout,
});

let N, X, Y, H;
let k = 0;
let precautionCount = 0;
const precautionLimit = 1e5;
let totalCount = 0;
let all = [];
let filtered = [];
let curr = "";
const nuc = "ACGT";
let reverseComplement = new Map();
let prohibited = new Map();

function initReverseComplement() {
  reverseComplement.set("A", "T");
  reverseComplement.set("T", "A");
  reverseComplement.set("G", "C");
  reverseComplement.set("C", "G");
}

function getReverseComplement(s) {
  let res = "";
  for (let a of s) {
    res += reverseComplement.get(a);
  }
  return res;
}

function generateRandomDNA(N) {
  const DNA = "ATCG";
  let randomDNA = "";
  for (let i = 0; i < N; ++i) {
    let randomIndex = Math.floor(Math.random() * DNA.length);
    randomDNA += DNA[randomIndex];
  }
  return randomDNA;
}

function costFunction(dna) {
  let uniqueCharacters = 0;
  for (let base of "ATCG") {
    if (dna.includes(base)) {
      uniqueCharacters++;
    }
  }
  return uniqueCharacters;
}

function stochasticLocalSearch(N, maxIterations) {
  let currentDNA = generateRandomDNA(N);

  for (let iteration = 1; iteration <= maxIterations; ++iteration) {
    let currentCost = costFunction(currentDNA);
    let neighborDNA = generateRandomDNA(N);
    let neighborCost = costFunction(neighborDNA);

    if (neighborCost < currentCost) {
      currentDNA = neighborDNA;
    }
  }

  return currentDNA;
}

function generateRandom() {
  for (let i = 0; i < 100000; i++) {
    all.push(stochasticLocalSearch(N, Math.floor(Math.random() * 15)));
  }
}

function generateCombinations() {
  if (curr.length === N) {
    all.push(curr);
    return;
  }

  if (totalCount >= 1e5) return;

  let st = Math.floor(Math.random() * 4);
  for (let i = 0; i < 4; i++) {
    let n = nuc[(i + st) % 4];
    curr += n;
    generateCombinations();
    curr = curr.slice(0, -1);
  }
}

function check(degree, s) {
  let n = s.length;
  for (let i = 0; i < n; i++) {
    let startOne = i;
    let startTwo = degree + i;
    let temp1 = s.substring(startOne, startOne + degree);
    let temp2 = "";
    if (startTwo < n) {
      temp2 = s.substring(startTwo, startTwo + degree);
    }
    if (temp1 === temp2) return false;
  }
  return true;
}

function homopolymerFree(s, limit) {
  for (let degree = 1; degree <= limit; degree++) {
    if (!check(degree, s)) return false;
  }
  return true;
}

function secondaryStructureFree(s, threshold) {
  let rc = s
    .split("")
    .reverse()
    .map((c) => reverseComplement.get(c))
    .join("");
  let store = N;
  N = s.length;
  let maxi = -1;
  let dp = Array.from({ length: N + 1 }, () => Array(N + 1).fill(0));

  for (let i = N; i >= 0; i--) {
    for (let j = N; j >= 0; j--) {
      if (i === N || j === N) {
        dp[i][j] = 0;
        continue;
      }
      if (s[i] === rc[j]) {
        dp[i][j] = 1 + dp[i + 1][j + 1];
        maxi = Math.max(maxi, dp[i][j]);
      } else {
        dp[i][j] = 0;
      }
    }
  }

  N = store;
  if (maxi > threshold) {
    return false;
  }
  return true;
}

function balancedGC(s) {
  let mp = new Map();
  for (let a of s) {
    mp.set(a, (mp.get(a) || 0) + 1);
  }
  return mp.get("A") + mp.get("T") === mp.get("G") + mp.get("C");
}

function findEditDis(s, t) {
  let n = s.length;
  let m = t.length;
  let dp = Array.from({ length: n + 1 }, () => Array(m + 1).fill(0));

  for (let i = 0; i <= n; i++) {
    for (let j = 0; j <= m; j++) {
      if (i === 0 || j === 0) {
        dp[i][j] = i + j;
      } else if (s[i - 1] === t[j - 1]) {
        dp[i][j] = dp[i - 1][j - 1];
      } else {
        dp[i][j] = Math.min(
          Math.min(
            1 + dp[i - 1][j], // Insert.
            1 + dp[i][j - 1]
          ), // Remove.
          1 + dp[i - 1][j - 1] // Replace.
        );
      }
    }
  }

  return dp[n][m];
}

let valid = false;

function getXCodewordsPrint(allcurr, i, curr) {
  if (valid === true) return;
  if (precautionCount >= precautionLimit) {
    return;
  }
  if (i === 0) {
    for (let a of curr) {
      for (let b of curr) {
        precautionCount++;
        if (a !== b && findEditDis(a, b) < H) return;
      }
    }
    valid = true;
    for (let a of curr) console.log(a);
    console.log("\n\n");
    return;
  }
  if (valid === true) return;
  for (let a of allcurr) {
    let flag = 0;
    for (let b of curr) {
      precautionCount++;
      if (findEditDis(a, b) < H) {
        flag = 1;
        break;
      }
    }
    if (flag) {
      continue;
    }
    curr.push(a);
    getXCodewordsPrint(allcurr, i - 1, curr);
    curr.pop();
  }
  return;
}

function getXCodewords(allcurr, i, curr) {
  if (valid === true) return;
  if (precautionCount >= precautionLimit) {
    return;
  }
  if (i === 0) {
    for (let a of curr) {
      for (let b of curr) {
        precautionCount++;
        if (a !== b && findEditDis(a, b) < H) return;
      }
    }
    valid = true;
    return;
  }
  if (valid === true) return;
  for (let a of allcurr) {
    let flag = 0;
    for (let b of curr) {
      precautionCount++;
      if (findEditDis(a, b) < H) {
        flag = 1;
        break;
      }
    }
    if (flag) {
      continue;
    }
    curr.push(a);
    getXCodewords(allcurr, i - 1, curr);
    curr.pop();
  }
  return;
}

rl.question("Enter the size of codeword (N): ", (size) => {
  N = parseInt(size);

  rl.question(
    "Enter X such that you want to prevent homopolymers up to X degrees: ",
    (x) => {
      X = parseInt(x);

      rl.question(
        "Enter the maximum length allowed for the patterns for secondary structures (Y): ",
        (y) => {
          Y = parseInt(y);

          rl.question(
            "Enter the minimum Edit / Edit distance required (H): ",
            (h) => {
              H = parseInt(h);

              initReverseComplement();
              totalCount = 0;

              if (N <= 10) generateCombinations();
              else generateRandom();

              for (let a of all) {
                if (homopolymerFree(a, X) && balancedGC(a)) {
                  if (secondaryStructureFree(a, Y)) {
                    prohibited.set(
                      getReverseComplement(a),
                      (prohibited.get(getReverseComplement(a)) || 0) + 1
                    );
                    if (prohibited.get(a) === 0) {
                      filtered.push(a);
                    }
                  }
                }
              }

              let siz = filtered.length;
              for (let i = 0; i < siz - 1; i++) {
                let j = i + Math.floor(Math.random() * (siz - i));
                [filtered[i], filtered[j]] = [filtered[j], filtered[i]];
              }

              let curr = [];
              let low = 1;
              let high = 300;

              while (low <= high) {
                let mid = Math.floor((low + high) / 2);
                console.log(`Checking ${mid}`);
                valid = false;
                precautionCount = 0;
                getXCodewords(filtered, mid + 1, curr);
                curr = [];
                if (valid) {
                  low = mid + 1;
                  console.log(` - Positive`);
                } else {
                  high = mid - 1;
                  console.log(` - Negative`);
                }
              }

              valid = false;
              let mid = Math.floor((low + high) / 2);
              precautionCount = 0;
              getXCodewordsPrint(filtered, mid, curr);

              if (!valid) {
                console.log(
                  "Sorry, no codewords can be formed with these requirements."
                );
              } else {
                console.log((low + high) / 2);
                console.log("You are having information density of:");
                let ans = Math.log2((low + high) / 2) / N;
                console.log(ans);
              }

              rl.close();
            }
          );
        }
      );
    }
  );
});
