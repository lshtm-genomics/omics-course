# Task

!!! Warning

    Be sure to use the `linux_hunt` Codespace when initiating VM for this practical.

## Introduction

Welcome to the **Linux Treasure Hunt**, a practical exercise designed to strengthen your command-line navigation and file manipulation skills. In real-world bioinformatics and data science work, much of your time is spent moving through complex directory structures, inspecting files, running scripts, and organising outputs efficiently. This task simulates that environment in a fun, exploratory way.

You will begin with a set of hidden clues inside a Linux directory. Each clue will guide you to the next location or action. By using core Linux commands, you will progress step-by-step towards uncovering the final treasure.

This exercise focuses on developing fluency with essential commands such as:
`ls`, `cd`, `cat`, `nano`, `mv`, `rm`, `mkdir`, and executing shell scripts.

---

## Task

A mysterious directory structure has been prepared for you in a Codespace environment. Somewhere within it lies a hidden treasure, but you must follow a sequence of clues to find it.

### Starting Point

Open the Linux test Codespace and begin here:

```bash
cd linux_hunt/
```
Your journey starts by reading the first clue:
```bash
cat clue1.txt
```
You will see a message similar to:
```bash
Welcome to the hunt, find the island and find the treasure. Remember to run bash scripts (.sh files) run bash "filename"
```
### Structure of the Hunt

The hunt is split into stages:

1. The Ship (starting area)
2. The Island
3. The Strange Tree
4. The Hollow Cave
5. The Hidden Chest

### Skills You Will Practice

- File navigation (`cd`, `ls`)
- File inspection (`cat`, `nano`)
- File manipulation (`mv`, `rm`, `touch`)
- Script execution (`bash`)
- Text searching (`grep`)
- Debugging broken workflows


### Challenge Rules

- Some files may regenerate after deletion
- Some clues may appear inside large log files
- Some scripts may fail until required conditions are met
- The same file may be required multiple times in different states

### End Goal

The hunt is complete when you successfully retrieve the final file:

treasure.txt

Inside it will contain the final message confirming completion of the challenge.


### Hints

If you are stuck:
- Re-read the last clue carefully
- Check file contents again using `cat`
- Search large files using `grep`
