// created by Richard Boyne rmb115@ic.ac.uk (last edit 20/02/19)
#pragma once
#include <map>
#include <iostream>
#include <string>
using namespace std;

// to use this menu create a class which inherits it which must have
	// - a constructor which sets the command options in command_map
	// - optionally set the same keys to a description in the descriptions
	// - optionally set the intro string printed at the start of the program
	// - there are a set of enter functions to help with standard use entries
class menu {
public:

	string intro = "No intro text given";  // text printed at the start
	string input;
	map<const string, int> command_map;  // maps command string to switch parameter
	map<const string, string> desctriptions;  // gives command descriptions


		// switch board of code to run on each command
	virtual void commands(int opt) = 0;


		// called this to beigin using the menu
	void start() {
			// repeatedly call get_input untill it returns 0 to quit
		while (true) {
			// print the intro texts
			cout << intro << '\n';
			cout << "=========================================================" << '\n'
				 << "                           Menu                          " << '\n'
				 << "=========================================================" << '\n'
				 << "\t Command \t Description \n" << endl;
			for (auto const& item : desctriptions)
				cout << '\t' << item.first << "\t" << item.second << '\n';
			cout << '\n'
				 << '\t' << "exit" << "\t\t" << "quit the program" << '\n'
				 << '\t' << "restart" << "\t\t" << "restart the program" << '\n'
				 << '\n' << "=========================================================" << "\n\n";


				// call UI which either returns 0 to quit or 1 to restart
			if (this->UI() == 0)
				return;
			else
				cout << "restarting \n"
				<< "--------------------------------------------------------- \n\n";
		}

	}


		// command to be overwritten - called before restarting
	virtual void clean() {
		return;
	}


		// user interface - repeatedly asks for commands
	int UI()
	{
			// repeatedly get inputs
		while (true)
		{
				// get the command
			cout << "Command: ";
			getline(cin, this->input);

				// if this command is specified in command_map
			if (command_map.count(input))
				commands(command_map[this->input]);

				// if user wants to quit
			else if (input == "exit") {
				this->clean();
				return 0;
			}

				// if user wants to restart
			else if (input == "restart") {
				this->clean();
				return 1;
			}

				// any other command is unknown
			else
				cout << "Command Unknown";

				// always start a new line
			cout << "\n\n";
		}
	}


		// simple interger entering function
	int enter_int(int min = INT32_MIN, int max = INT32_MAX) {
		string buffer;
		int n;
		while (true) {
				// put input in string buffer (as ' ' break cin)
			getline(cin, buffer);
			try {
				n = stoi(buffer);
				if (n < min || n > max) // check in range
					cout << "Error: Input must be within [" << min << ", " << max << "]\n";
				else
					return n;
			}
				// if could not convert to integer
			catch (const std::invalid_argument) {
				cout << "Error: Input not an integer\n";
			}
			cout << "\nPlease try again: ";
		}
	}


		// simple string entering function
	string enter_string(string must_contain = "", string cant_contain = "") {
		string buffer;
		while (true) {

			getline(cin, buffer);

			if (buffer.find(must_contain) == string::npos)
				cout << "Error: string must contain '" << must_contain << "'\n";
			else if (buffer.find(cant_contain) != string::npos)
				cout << "Error: string can't have '" << cant_contain << "'\n";
			else
				return buffer;
			cout << "\nPlease try again: ";
		}
	}


		// simple enter from an array of options
	template<class T>
	T enter_options(T* ptr, int length) {
		
		T buffer;
		bool found = 0;

		while (true) {

			cin >> buffer;
			cin.ignore(100, '\n');  // igore everything else till the new line
			
				// check input worked
			if (!cin) {
				cout << "Error: input failed to cast into vairable";
				continue;
			}

				// check if input is an option
			for (int i = 0; i < length; i++)
				if (ptr[i] == buffer)
					found = 1;
			if (found)
				return buffer;

			cout << "Error: input must be one of: ";
			for (int i = 0; i < length; i++)
				cout << ptr[i] << ", ";

			cout << "\n\nPlease try again: ";
		}
	}

		
		// simple enter bool option
	bool enter_bool() {
		
		char options[] = { 'y', 'n' };
		char entry = this->enter_options<char>(options, 2);

		if (entry == 'y')
			return 1;
		else
			return 0;
	}

};