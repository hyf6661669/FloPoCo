#ifndef ExpressionParserData_H
#define ExpressionParserData_H

#include <string>
#include <vector>
#include <map>


namespace flopoco{


class DFGNode;

/**
 * @brief The Node class
 * the main explanation needed is about the signification of the
 */

enum NodeTypeEnum {EntryVar,OpAssign,OpConst,OpAdd,OpSub,OpMult,OpDiv,OpLog,OpExp};


class DFGNode{

public:
    DFGNode(){
        depth=-1;
    }


    NodeTypeEnum type_id;//id du type de node
    int depth;

    std::vector<DFGNode*> parents;
    std::vector<DFGNode*> childrens;

    int nbr_childrens(){
        return (int)childrens.size();
    }

    int nbr_parents(){
        return (int)parents.size();
    }

    void rmv_child(DFGNode* node){
        for (int i=0;i<(int)childrens.size();++i){
            if (childrens[i]==node){
                childrens.erase(childrens.begin()+i);
                --i;
            }
        }
    }

    void rmv_parent(DFGNode* node){
        for (int i=0;i<(int)parents.size();++i){
            if (parents[i]==node){
                parents.erase(parents.begin()+i);
                --i;
            }
        }
    }


};



struct Program{
    std::vector<DFGNode*> return_heads;
    std::map<std::string,DFGNode*> variables;

    void propagate_depth(DFGNode* node, int depth){
        if (node->depth<depth){
            node->depth=depth;
            std::cout<<node->type_id<<"::"<<depth<<std::endl;
        }else{
            std::cout<<"depth already achieved"<<std::endl;
            return;
        }

        for (int i=0;i<(int)node->childrens.size();++i){
            propagate_depth(node->childrens[i],depth+1);
        }
    }

    void init_depth(){
        for (int i=0;i<(int)return_heads.size();++i){
            std::cout<<"starting new tree"<<std::endl;
            propagate_depth(return_heads[i],0);
        }
    }
};



//this is how to read a Number*
/*
  Number* nbr;
  nbr= new Integer();

  Integer* it= dynamic_cast<Integer*>(nbr);
  //possible to acces the integer members
  if (it!=NULL){
    use the it variable
  }

  FloatingPointNumber* fp=dynamic_cast<FloatingPointNumber*>(nbr);
  //impossible to acces the FloatingPointNumber members the fp variable was set to NULL
  if (fp!=NULL){
    won't be here
  }
 */


class Number{

};

class FloatingPointNumber:public Number{
public:
    int mantissa_size;
    long mantissa;

    int exponent_size;
    long exponent;
};

class FixedPointNumber:public Number{
    ///TODO implement FixedPointNumbers
};

class Integer:public Number{
public:
    bool is_signed;
    int size;
    long value;
};


class EntryVariable:public DFGNode{
public:
    Number* nbr;

    EntryVariable():DFGNode(){
        type_id=EntryVar;
    }
};

class OperatorConstant: public DFGNode{
public:
    NodeTypeEnum operator_type;//indicate *  or +  or /  or -  if it's just a constant number this value indicate (OpConst)

    Number* nbr;

    OperatorConstant():DFGNode(){
        type_id=OpConst;
    }

};

}

#endif
