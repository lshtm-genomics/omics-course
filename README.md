# omics-course

## Contribution

Here are instructions on how to update/add new tutorials

### Setting up

It is best to do this in a fresh conda environment which you can make with

```
mamba create -n mkdocs mkdocs-material
```

Then on your computer clone the repository:

```
git clone https://github.com/lshtm-genomics/omics-course
```

Checkout the `dev` branch:

```
git checkout dev
```

Then go into the folder and run mkdocs serve

```
mkdocs serve
```

### Adding new tutorials

If you need add a new tutorial just create a new markdown file in `docs/`. If you are running mkdocs serve your changes will automatically force the browser to refresh.

If you make some changes then just use:

```
git add docs/new-tutorial.md
git commit -am 'added new tutorial'
git push
```

#### Images

All images can be placed using the following notation

```
![mapping_1](../img/Mapping_1.jpg)
```

#### Questions

If there is a question that you want the participant to think about you can format it like this:

```
!!! question "Question 1"
    Is this a question?
```

If your question has a specific answer you can use the following formatting

```
!!! question
    === "Question 1"
        Is this a question?
    === "Answer 1"
        This is an answer
```

#### Information

Tips and information can be inserted using the following:

```
!!! info
    This is some information
```

#### Terminal output

If you want to put in some expected output from the terminal that the participants can compare to then this will work:

```
!!! terminal "Terminal output"
    ```
    This is osme terminal output
    ```
```
#### Code

Code can be inserted with the the triple backticks (without the backslash):

```
\```
some_command -a parameter1
\```
```


