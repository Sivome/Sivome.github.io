---
title: "import pdb; pdb.set_trace() <EOB>"
date: '2019-03-18'
layout: post
categories: Proteomics
---

Take my word for this. This will help you for years to come if you are in data science.
By the way,
```python
import pdb; pdb.set_trace()
```

Just include it anywhere in your python script for debugging.  It will put the script in the debug mode. In the debug mode, other than the general features of viewing the current variables etc, I always add new variables, come up with new functionality that I can later implement in the original script after quitting the model.  

It helped me a ton in the last couple of scripts I wrote, and I thought **it** deserves a blog post by itself. My original intention is to just include the command in the title, and so I ended it with EOB (End of Blog).  Do NOT include <EOB> in the debug mode :-)

More info here: https://docs.python.org/2/library/pdb.html  
